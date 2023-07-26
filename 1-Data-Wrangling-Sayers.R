
# script designed to join various Hg databases and match species to appropriate
# functional traits

#---------------------- LOADING/MERGING THE DATA -------------------------------
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(janitor)

# resolve namespace conflicts and creating necessary functions
select <- dplyr::select
"%nin%" <- Negate("%in%")

# Cumulative TRACE database = BRI + SDZWA + Duke + CINCIA + Shrum data
TRACEData <- read.csv("TRACE_Database_26Jul2023.csv", na.strings = c("",".","NA")) %>% 
  # removing captive birds from data set
  filter(Site_Name %nin% c("Belize Zoo", "Belize Raptor Center")) %>% 
  # following the assumption that 95% of THg in feathers is MeHg, we can effectively compare
  # MeHg concentrations with THg concentrations
  mutate(Flank_MeHg_ppm = Flank_MeHg_ppb / 1000) %>% # converting from ppb to ppm
  # joining flank THg and MeHg concentrations together
  unite("Flank_Hg_ppm", c(Flank_Hg_ppm, Flank_MeHg_ppm), na.rm = TRUE, remove = F, sep = "") %>%
  transform(Flank_Hg_ppm = as.numeric(Flank_Hg_ppm)) %>% # converting to numeric for later computation
  # joining breast, flank, and back feather concentrations together to create a body feather category
  unite("Body_Hg_ppm", c(Breast_Hg_ppm, Flank_Hg_ppm, Back_Hg_ppm), na.rm = TRUE, remove = F, sep = "") %>% 
  transform(Body_Hg_ppm = as.numeric(Body_Hg_ppm)) %>% # converting to numeric for later computation
  # excluding samples below the lower detection limit of the Hg analyzer
  mutate(Blood_Hg_ppm = ifelse(Blood_Hg_ppm < 0.001, NA, Blood_Hg_ppm),
         Tail_Hg_ppm = ifelse(Tail_Hg_ppm < 0.001, NA, Tail_Hg_ppm),
         Breast_Hg_ppm = ifelse(Breast_Hg_ppm < 0.001, NA, Breast_Hg_ppm),
         Flank_Hg_ppm = ifelse(Flank_Hg_ppm < 0.001, NA, Flank_Hg_ppm),
         Back_Hg_ppm = ifelse(Back_Hg_ppm < 0.001, NA, Back_Hg_ppm),
         Body_Hg_ppm = ifelse(Body_Hg_ppm < 0.001, NA, Body_Hg_ppm)) %>%
  # excluding sparse historical feather samples from museums
  filter(Year > 2006) %>%  
  # creating a seasonal column distinguishing wet and dry season by region
  # Belize wet season: June through December
  mutate(Season = if_else(Country == "Belize" & Month %in% c(6:12), "Wet",
                          # western Mexico wet season: June through October
                          if_else(Country == "Mexico" & Month %in% c(6:10), "Wet",
                                  # Puerto Rico wet season: May through December
                                  if_else(Country == "Puerto Rico" & Month %in% c(4:12), "Wet",
                                          # Dominican Republic wet season: May through December
                                          if_else(Country == "Dominican Republic" & Month %in% c(5:11), "Wet",
                                                  # Nicaragua wet season: May through October                    
                                                  if_else(Country == "Nicaragua" & Month %in% c(5:11), "Wet",
                                                          # Costa Rica wet season: May through October       
                                                          if_else(Country == "Costa Rica" & Month %in% c(5:11), "Wet",
                                                                  # Panama wet season: May through November      
                                                                  if_else(Country == "Panama" & Month %in% c(5:11), "Wet", 
                                                                          # Madre de Dios, Peru wet season: October through April      
                                                                          if_else(Country == "Peru" & Month %in% c(1:4, 10:12), "Wet",
                                                                                  # treat all other seasons as "dry"
                                                                                  "Dry")))))))))


##### ADDING ORDER AND FAMILY USING eBIRD/CLEMENTS v2019 CLASSIFICATIONS #####

# not sure why the .csv file isn't available yet
taxa <- read.csv("NEW_eBird-Clements-v2022-integrated-checklist-October-2022.csv") %>%
  rename(Species_Latin_Name = SCI_NAME, Order = ORDER1, Family = FAMILY) %>% # renaming key column names
  select(Species_Latin_Name, Order, Family)

# Which species were not joined due to taxonomic changes?
unjoined <- anti_join(TRACEData, taxa, by = "Species_Latin_Name")
unique(unjoined$Species_Latin_Name) # This should be 0 -- It is!


##### CLASSIFYING SPECIES VIA HABITAT USING PARKER et al. (1996) CRITERIA #####

adata <- read.csv("Parker_Stotz_Fitzpatrick_1996/adata.csv") %>% # neotropical breeder traits
  mutate(Migratory_Status = "Resident")
cdata <- read.csv("Parker_Stotz_Fitzpatrick_1996/cdata.csv") %>%  # neotropical migrant nonbreeding traits
  mutate(Migratory_Status = "Full migrant")
ddata <- read.csv("Parker_Stotz_Fitzpatrick_1996/ddata.csv") %>% # neotropical migrant breeding traits
  mutate(Migratory_Status = "Partial migrant")
edata <- read.csv("Parker_Stotz_Fitzpatrick_1996/edata.csv") %>% # austral migrant traits
  mutate(Migratory_Status = "Austral migrant")

# combing relevant data frames for our purposes (a, c, and e) by rows without duplicating column headings
parker <- full_join(cdata, ddata) %>%
  full_join(edata) %>%
  full_join(adata) %>%
  unite(Species_Latin_Name, c("GENUS", "SPECIES"), sep = " ", remove = F) %>% 
  # this function determines which neotropical "residents" are duplicated as partial
  # migrants in "ddata" or "austral migrants in "edata", it then returns only the
  #data associated with the resident entry (e.g. Egretta alba -- old name for Great Egret)
  distinct(Species_Latin_Name, .keep_all = T) %>% 
  # fixing migratory classes of certain species
  mutate(Migratory_Status = if_else(Species_Latin_Name == "Elanoides forficatus", "Partial migrant", Migratory_Status)) %>% # Swallow-tailed Kite
  mutate(Migratory_Status = if_else(Species_Latin_Name == "Haematopus palliatus", "Partial migrant", Migratory_Status)) %>% # American Oystercatcher
  mutate(Migratory_Status = if_else(Species_Latin_Name == "Himantopus mexicanus", "Partial migrant", Migratory_Status)) %>% # Black-necked Stilt
  mutate(Migratory_Status = if_else(Species_Latin_Name == "Mycteria americana", "Partial migrant", Migratory_Status)) %>% # Wood Stork
  mutate(Migratory_Status = if_else(Species_Latin_Name == "Progne chalybea", "Partial migrant", Migratory_Status)) %>% # Gray-breasted Martin
  mutate(Migratory_Status = if_else(Species_Latin_Name == "Butorides virescens", "Partial migrant", Migratory_Status)) # Green Heron

# Which species were not joined due to taxonomic changes or subspecies?
unjoined <- anti_join(TRACEData, parker, by = "Species_Latin_Name")
unique(unjoined$Species_Latin_Name) # 123 taxa had issues

# Manually changing Parker et al. dataset to accommodate for changed species names, lumps,
# and splits since 1996 using Birds of the World taxonomic criteria. Only changing species
# relevant to this Hg study (for right now at least). Jacob Socolar made excellent headway
# on this, but some of his changes are already outdated given recent taxonomic updates.
# After pitching the project to the BOW team, they are currently discussing creating
# a resource similar to the Parker criteria but for all bird species — so, this work
# should at least serve in the interim of that endeavor being completed.

# NOTE: some species may appear twice with different life history information in the
# Parker file (ie. Yellow Warbler) because subspecies may be resident and migratory

# Simple name changes since 1996 -- this does NOT fix the order and family changes since 1996, 
# but we will handle that with the eBird/Clements database
# Old species name in Parker <- new species name in BOW

# New name is being assigned to old Parker name
parker <- parker %>%
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Actitis macularia", "Actitis macularius", Species_Latin_Name)) %>% # Spotted Sandpiper
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Cacicus holosericeus", "Amblycercus holosericeus", Species_Latin_Name)) %>% # Yellow-billed Cacique
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Atlapetes brunneinucha", "Arremon brunneinucha", Species_Latin_Name)) %>% # Chestnut-capped Brush-Finch
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Arremon (taciturnus) taciturnus", "Arremon taciturnus", Species_Latin_Name)) %>% # Pectoral Sparrow
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Ardeola ibis", "Bubulcus ibis", Species_Latin_Name)) %>% # Cattle Egret
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Buteo magnirostris", "Rupornis magnirostris", Species_Latin_Name)) %>% # Roadside Hawk
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Butorides (striatus) virescens", "Butorides virescens", Species_Latin_Name)) %>% # Green Heron
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Wilsonia pusilla", "Cardellina pusilla", Species_Latin_Name)) %>% # Wilson's Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Ceryle torquata", "Megaceryle torquata", Species_Latin_Name)) %>% # Ringed Kingfisher
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Chlorospingus ophthalmicus", "Chlorospingus flavopectus", Species_Latin_Name)) %>% # Common Chlorospingus
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Contopus borealis", "Contopus cooperi", Species_Latin_Name)) %>% # Olive-sided Flycatcher
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Piaya minuta", "Coccycua minuta", Species_Latin_Name)) %>% # Little Cuckoo
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Scardafella inca", "Columbina inca", Species_Latin_Name)) %>% # Inca Dove
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Columbina (talpacoti) talpacoti", "Columbina talpacoti", Species_Latin_Name)) %>% # Ruddy Ground Dove
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Passerina parellina", "Cyanocompsa parellina", Species_Latin_Name)) %>% # Blue Bunting
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Dendrocincla (fuliginosa) fuliginosa", "Dendrocincla fuliginosa", Species_Latin_Name)) %>% # Plain-brown Woodcreeper
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Xiphorhynchus picus", "Dendroplex picus", Species_Latin_Name)) %>% # Straight-billed Woodcreeper
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Pipra pipra", "Dixiphia pipra", Species_Latin_Name)) %>% # White-crowned Manakin
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Geothlypis formosus", "Geothlypis formosa", Species_Latin_Name)) %>% # Kentucky Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Glaucis hirsuta", "Glaucis hirsutus", Species_Latin_Name)) %>% # Rufous-breasted Hummingbird
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Myrmeciza fortis", "Hafferia fortis", Species_Latin_Name)) %>% # Sooty Antbird
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Helmitheros vermivorus", "Helmitheros vermivorum", Species_Latin_Name)) %>% # Worm-eating Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Icterus (galbula) galbula", "Icterus galbula", Species_Latin_Name)) %>% # Baltimore Oriole
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Icterus (spurius) spurius", "Icterus spurius", Species_Latin_Name)) %>% # Orchard Oriole
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Myrmotherula hauxwelli", "Isleria hauxwelli", Species_Latin_Name)) %>% # Plain-throated Antwren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Vermivora peregrina", "Leiothlypis peregrina", Species_Latin_Name)) %>% # Tennessee Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Leptotila (rufaxilla) plumbeiceps", "Leptotila plumbeiceps", Species_Latin_Name)) %>% # Gray-headed Dove
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Leptotila (rufaxilla) rufaxilla", "Leptotila rufaxilla", Species_Latin_Name)) %>% # Gray-fronted Dove
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Larus atricilla", "Leucophaeus atricilla", Species_Latin_Name)) %>% # Laughing Gull
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Myrmotherula fulviventris", "Epinecrophylla fulviventris", Species_Latin_Name)) %>% # Checker-throated Stipplethroat
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Gymnopithys salvini", "Oneillornis salvini", Species_Latin_Name)) %>% # White-throated Antbird
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Onychorhynchus (coronatus) coronatus", "Onychorhynchus coronatus", Species_Latin_Name)) %>% # Royal Flycatcher
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Geothlypis agilis", "Oporornis agilis", Species_Latin_Name)) %>% # Connecticut Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Vermivora gutturalis", "Oreothlypis gutturalis", Species_Latin_Name)) %>% # Flame-throated Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Oryzoborus (angolensis) funereus", "Sporophila funerea", Species_Latin_Name)) %>% # Thick-billed Seed-Finch
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Seiurus motacilla", "Parkesia motacilla", Species_Latin_Name)) %>% # Lousiana Waterthrush
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Seiurus noveboracensis", "Parkesia noveboracensis", Species_Latin_Name)) %>% # Northern Waterthrush
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Thryothorus (maculipectus) maculipectus", "Pheugopedius maculipectus", Species_Latin_Name)) %>% # Spot-breasted Wren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Pipra mentalis", "Ceratopipra mentalis", Species_Latin_Name)) %>% # Red-capped Manakin
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Platyrinchus (mystaceus) cancrominus", "Platyrinchus cancrominus", Species_Latin_Name)) %>% # Stub-tailed Spadebill
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Rhynchocyclus (brevirostris) brevirostris", "Rhynchocyclus brevirostris", Species_Latin_Name)) %>% # Eye-ringed Flatbill
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Seiurus aurocapillus", "Seiurus aurocapilla", Species_Latin_Name)) %>% # Ovenbird
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Parula americana", "Setophaga americana", Species_Latin_Name)) %>% # Northern Parula
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Wilsonia citrina", "Setophaga citrina", Species_Latin_Name)) %>% # Hooded Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Dendroica discolor", "Setophaga discolor", Species_Latin_Name)) %>% # Prairie Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Dendroica magnolia", "Setophaga magnolia", Species_Latin_Name)) %>% # Magnolia Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Dendroica pensylvanica", "Setophaga pensylvanica", Species_Latin_Name)) %>% # Chestnut-sided Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Dendroica (petechia) aestiva", "Setophaga petechia", Species_Latin_Name)) %>% # Yellow Warbler
  # NOTE: because BRI did not collect subspecies data, I cannot definitely determine the subspecies
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Dendroica tigrina", "Setophaga tigrina", Species_Latin_Name)) %>% # Cape May Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Oryzoborus (angolensis) angolensis", "Sporophila angolensis", Species_Latin_Name)) %>% # Chestnut-bellied Seed-Finch
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Sporophila (torqueola) morelleti", "Sporophila morelleti", Species_Latin_Name)) %>% # Morelet's Seedeater
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Tangara larvata", "Stilpnia larvata", Species_Latin_Name)) %>% # Golden-hooded Tanager
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Synallaxis (gujanensis) gujanensis", "Synallaxis gujanensis", Species_Latin_Name)) %>% # Plain-crowned Spinetail
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Thamnophilus (doliatus) doliatus", "Thamnophilus doliatus", Species_Latin_Name)) %>% # Barred Antshrike
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Thryothorus rufalbus", "Thryophilus rufalbus", Species_Latin_Name)) %>% # Rufous-and-white Wren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Thryothorus semibadius", "Cantorchilus semibadius", Species_Latin_Name)) %>% # Riverside Wren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Troglodytes (aedon) aedon", "Troglodytes aedon", Species_Latin_Name)) %>% # House Wren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Vermivora pinus", "Vermivora cyanoptera", Species_Latin_Name)) %>% # Blue-winged Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Vireo (olivaceus) olivaceus", "Vireo olivaceus", Species_Latin_Name)) %>% # Red-eyed Vireo
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Myrmeciza hemimelaena", "Sciaphylax hemimelaena", Species_Latin_Name)) %>% # Chestnut-tailed Antbird (this  species was split but then lumped back)
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Myrmeciza hyperythra", "Myrmelastes hyperythrus", Species_Latin_Name)) %>% # Plumbeous Antbird
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Oryzoborus (maximiliani) maximiliani", "Sporophila atrirostris", Species_Latin_Name)) %>% # Black-billed Seed-Finch
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Thryothorus (genibarbis) genibarbis", "Pheugopedius genibarbis", Species_Latin_Name)) %>% # Moustached Wren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Veniliornis passerinus", "Dryobates passerinus", Species_Latin_Name)) %>%  # Little Woodpecker
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Xiphorhynchus (spixii) elegans", "Xiphorhynchus elegans", Species_Latin_Name)) %>% # Elegant Woodcreeper
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Momotus (momota) momota", "Momotus momota", Species_Latin_Name)) %>% # Amazonian Motmot
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Myrmeciza melanoceps", "Akletos melanoceps", Species_Latin_Name)) %>% # White-shouldered Antbird
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Pipra coronota", "Lepidothrix coronata", Species_Latin_Name)) %>% # Blue-crowned Manakin
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Corythopis torquata", "Corythopis torquatus", Species_Latin_Name)) %>% # Ringed Antpipit
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Thryothorus leucotis", "Cantorchilus leucotis", Species_Latin_Name)) %>% # Buff-breasted Wren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Donacobius atricapillus", "Donacobius atricapilla", Species_Latin_Name)) %>% # Black-capped Donacobius
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Myrmeciza atrothorax", "Myrmophylax atrothorax", Species_Latin_Name)) %>% # Black-throated Antbird
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Otus watsonii", "Megascops watsonii", Species_Latin_Name)) %>% # Tawny-bellied Screech-Owl
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Otus choliba", "Megascops choliba", Species_Latin_Name)) %>% # Tropical Screech-Owl
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Pipra chloromeros", "Ceratopipra chloromeros", Species_Latin_Name)) %>% # Round-tailed Manakin
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Amazilia candida", "Chlorestes candida", Species_Latin_Name)) %>% # White-bellied Emerald
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Chlorostilbon canivetii", "Cynanthus canivetii", Species_Latin_Name)) %>% # Canivet's Emerald
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Glaucidium brasilianum", "Nannopterum brasilianum", Species_Latin_Name)) %>% # Neotropical Cormorant
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Hylocharis cyanus", "Chlorestes cyanus", Species_Latin_Name)) %>% # White-chinned Sapphire
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Pteroglossus beauharnaesii", "Pteroglossus beauharnaisii", Species_Latin_Name)) %>% # Curl-crested Aracari
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Tachyphonus luctuosus", "Loriotus luctuosus", Species_Latin_Name)) %>% # White-shouldered Tanager
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Amazilia lactea", "Chionomesa lactea", Species_Latin_Name)) %>% # Sapphire-spangled Emerald
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Veniliornis fumigatus", "Dryobates fumigatus", Species_Latin_Name)) %>% # Smoky-brown Woodpecker
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Icterus (dominicensis) prosthemelas", "Icterus prosthemelas", Species_Latin_Name)) %>% # Black-cowled Oriole
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Hylophilus decurtatus", "Pachysylvia decurtata", Species_Latin_Name)) %>% # Lesser Greenlet
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Leucopternis schistacea", "Buteogallus schistaceus", Species_Latin_Name)) %>% # Slate-colored Hawk
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Speotyto cunicularia", "Athene cunicularia", Species_Latin_Name)) # Burrowing Owl

# Which species were STILL not joined due to taxonomic changes or subspecies? — eventually will be 0
unjoined <- anti_join(TRACEData, parker, by = "Species_Latin_Name")
unique(unjoined$Species_Latin_Name)
# We have 35 taxa that have been split/lumped since 1996 or were just excluded from the df

# Manually adding species to accommodate for splits/lumps information using BOW criteria
# I am choosing to leave the original classifications in the df so that we can reference it later if necessary
# If I was unsuccessful in finding information on BOW, I carried that information down from the old species in Parker
# ^^ This is common for poorly studied neotropical species

parker <- parker %>% 
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "TROGLODYTIDAE",
          Species_Latin_Name = "Cantorchilus modestus", # Cabanis's Wren, formerly Thryothorus modestus
          GENUS = "Cantorchilus",
          SPECIES = "modestus",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "L", # Sensitivity
          STRAT = "U/M", # Foraging Strata (Chris added M)
          CNTAB = "LT", # Center of Abundance    
          REL = "C", # Relative Abundance
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 2000, # Maximum Elevation (Chris changed to 2000 m)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA, # Microhabitats
          NHAB = 3, # Number of Habitats
          NZOO = 4, # Number of Zoogeographic Regions
          HAB1 = "F7E", HAB2 = "F1E", HAB3 = "N14", HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "E", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = "E", F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = NA,
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = "Y",
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = "Y", CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "TROGLODYTIDAE",
          Species_Latin_Name = "Cantorchilus thoracicus", # Stripe-breasted Wren, formerly Thryothorus thoracicus (split into Cantorchilus leucopogon)
          GENUS = "Cantorchilus",
          SPECIES = "leucopogon",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "M", # Sensitivity
          STRAT = "U/M", # Foraging Strata
          CNTAB = "LT", # Center of Abundance    
          REL = "U/F", # Relative Abundance
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1100, # Maximum Elevation (Chris changed to 1100)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = "o", # Subregions
          MICRO = "S", # Microhabitats (Chris added S)
          NHAB = 3, # Number of Habitats
          NZOO = 3, # Number of Zoogeographic Regions (Chris added F3)
          HAB1 = "F1E", HAB2 = "F3", HAB3 = "F15", HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "E", F2 = NA, F3 = "Y", F4 = NA, F5 = NA, F6 = NA, F7 = NA, F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "Y",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = "Y",
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = "Y", MAH = "Y", GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "THAMNOPHILIDAE",
          Species_Latin_Name = "Cercomacroides tyrannina", # Dusky Antbird, formerly Cercomacra tyrannina
          GENUS = "Cercomacroides",
          SPECIES = "tyrannina",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "L", # Sensitivity
          STRAT = "U", # Foraging Strata
          CNTAB = "LT", # Center of Abundance    
          REL = "F/C", # Relative Abundance (Chris added F)
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1900, # Maximum Elevation (Chris changed to 1900)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = "V/A", # Microhabitats (Chris added V/A)
          NHAB = 2, # Number of Habitats
          NZOO = 7, # Number of Zoogeographic Regions (Chris changed)
          HAB1 = "F1E", HAB2 = "F15", HAB3 = NA, HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "E", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = NA, F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "Y",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = "Y", GCS = "Y", CDH = "Y", CHO = "Y",    
          EPC = NA, STP = NA, NAN = "Y", CAN = NA, SAN = NA, NSA = "Y", TEP = NA, AMN = "Y", AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "CARDINALIDAE",
          Species_Latin_Name = "Cyanoloxia cyanoides", # Blue-Black Grosbeak, formerly Passerina cyanoides
          GENUS = "Cyanoloxia",
          SPECIES = "cyanoides",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "M", # Sensitivity
          STRAT = "U", # Foraging Strata
          CNTAB = "LT", # Center of Abundance    
          REL = "C", # Relative Abundance (Chris changed to C from F)
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1250, # Maximum Elevation
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = "T", # Microhabitats (Chris added T)
          NHAB = 2, # Number of Habitats
          NZOO = 4, # Number of Zoogeographic Regions (Chris changed)
          HAB1 = "F1", HAB2 = "F15", HAB3 = NA, HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "Y", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = NA, F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "Y",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = "Y", CDH = NA, CHO = "Y",    
          EPC = "Y", STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = "Y", TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "CARDINALIDAE",
          Species_Latin_Name = "Cyanoloxia rothschildii", # Amazonian Grosbeak, formerly Passerina cyanoides
          GENUS = "Cyanoloxia",
          SPECIES = "rothschildii",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "M", # Sensitivity
          STRAT = "U", # Foraging Strata
          CNTAB = "LT", # Center of Abundance    
          REL = "F", # Relative Abundance
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1400, # Maximum Elevation (Chris changed to 1400)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA, # Microhabitats
          NHAB = 2, # Number of Habitats
          NZOO = 3, # Number of Zoogeographic Regions (Chris changed)
          HAB1 = "F1", HAB2 = "F15", HAB3 = NA, HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "Y", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = NA, F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "Y",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = "Y", TEP = NA, AMN = "Y", AMS = "Y",   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "FURNARIIDAE",
          Species_Latin_Name = "Dendrocolaptes sanctithomae", # Northern Barred-Woodcreeper, formerly Dendrocolaptes (certhia) certhia
          GENUS = "Dendrocolaptes",
          SPECIES = "sanctithomae",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "H", # Sensitivity
          STRAT = "U/M", # Foraging Strata
          CNTAB = "LT", # Center of Abundance    
          REL = "F", # Relative Abundance
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1800, # Maximum Elevation (Chris changed to 1800)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = "A", # Microhabitats
          NHAB = 7, # Number of Habitats (Chris added 6 HABs)
          NZOO = 7, # Number of Zoogeographic Regions (Chris removed AMN, AMS, ATL)
          HAB1 = "F1", HAB2 = "F15", HAB3 = "F14", HAB4 = "F8", HAB5 = "F7", HAB6 = "F4", HAB7 = "F11", #Habitat
          F1 = "Y", F2 = NA, F3 = NA, F4 = "Y", F5 = NA, F6 = NA, F7 = "Y", F8 = "Y",
          F9 = NA, F10 = NA, F11 = "Y", F12 = NA, F13 = NA, F14 = "Y", F15 = "Y",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = "Y", GCS = "Y", CDH = NA, CHO = "Y",    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = "Y", TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "THAMNOPHILIDAE",
          Species_Latin_Name = "Epinecrophylla amazonica", # Rio Madeira Stipplethroat, formerly known as Myrmotherula (haematonota) haematonota
          GENUS = "Epinecrophylla",
          SPECIES = "amazonica",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "H", # Sensitivity
          STRAT = "U", # Foraging Strata
          CNTAB = "LT", # Center of Abundance    
          REL = "F", # Relative Abundance
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 400, # Maximum Elevation (Chris changed to 400)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 2, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA, # Microhabitats
          NHAB = 1, # Number of Habitats
          NZOO = 1, # Number of Zoogeographic Regions
          HAB1 = "F1", HAB2 = NA, HAB3 = NA, HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "Y", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = NA, F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = NA,
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = NA, AMS = "Y",   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "FURNARIIDAE",
          Species_Latin_Name = "Furnarius leucopus", # Pale-legged Hornero, formerly Furnarius (leucopus) leucopus and Furnarius (leucopus) cinnamomeus
          GENUS = "Furnarius",
          SPECIES = "leucopus",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "L", # Sensitivity
          STRAT = "T", # Foraging Strata
          CNTAB = "LT", # Center of Abundance    
          REL = "U/C", # Relative Abundance (Chris added U)
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 2700, # Maximum Elevation (Chris changed to 400)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = "MAN", # Subregions
          MICRO = NA, # Microhabitats
          NHAB = 5, # Number of Habitats
          NZOO = 5, # Number of Zoogeographic Regions (Chris added EPC, N13, F15E, removed F2)
          HAB1 = "N14", HAB2 = "N13", HAB3 = "F15E", HAB4 = "F3", HAB5 = "F8", HAB6 = NA, HAB7 = NA, #Habitat
          F1 = NA, F2 = NA, F3 = "Y", F4 = NA, F5 = NA, F6 = NA, F7 = NA, F8 = "Y",
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "E",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = "Y", N14 = "Y",
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = "Y", STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = "Y", TEP = NA, AMN = "Y", AMS = "Y",   
          CSA = "Y", ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "THAMNOPHILIDAE",
          Species_Latin_Name = "Hypocnemis peruviana", # Peruvian-Warbling Antbird, formerly Hypocnemis cantator
          GENUS = "Hypocnemis",
          SPECIES = "peruviana",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "M", # Sensitivity
          STRAT = "U/M", # Foraging Strata
          CNTAB = "LT", # Center of Abundance    
          REL = "F/C", # Relative Abundance (Chris added F)
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1100, # Maximum Elevation
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = "A/T", # Microhabitats (Chris added A/T)
          NHAB = 3, # Number of Habitats
          NZOO = 2, # Number of Zoogeographic Regions
          HAB1 = "F1", HAB2 = "F2", HAB3 = "F15", HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "Y", F2 = "Y", F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = NA, F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "Y",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = "Y", AMS = "Y",   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "CORACIIFORMES",
          FAMILY = "MOMOTIDAE",
          Species_Latin_Name = "Momotus lessonii", # Lesson's Motmot, formerly Momotus (momota) momota
          GENUS = "Momotus",
          SPECIES = "lessonii",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "M", # Sensitivity
          STRAT = "U/M", # Foraging Strata
          CNTAB = "LT", # Center of Abundance    
          REL = "C", # Relative Abundance
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 2150, # Maximum Elevation (Changed from 1300)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA, # Microhabitats (Chris added A/T)
          NHAB = 6, # Number of Habitats
          NZOO = 1, # Number of Zoogeographic Regions
          HAB1 = "F1", HAB2 = "F4", HAB3 = "F15", HAB4 = "F8", HAB5 = "F7", HAB6 = "F2", HAB7 = NA, #Habitat
          F1 = "Y", F2 = "Y", F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = "Y", F8 = "Y",
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "Y",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = "Y", CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "THAMNOPHILIDAE",
          Species_Latin_Name = "Myrmelastes leucostigma", # Spot-winged Antbird, formely Percnostola leucostigma
          GENUS = "Myrmelastes",
          SPECIES = "leucostigma",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "H", # Sensitivity
          STRAT = "T/U", # Foraging Strata
          CNTAB = "HT", # Center of Abundance    
          REL = "C", # Relative Abundance
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1650, # Maximum Elevation (Changed from 1100)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = "S", # Microhabitats
          NHAB = 3, # Number of Habitats
          NZOO = 5, # Number of Zoogeographic Regions (Chris added F2)
          HAB1 = "F1", HAB2 = "F4", HAB3 = "F2", HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "Y", F2 = "Y", F3 = NA, F4 = "Y", F5 = NA, F6 = NA, F7 = NA, F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = NA,
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = "Y", CAN = "Y", SAN = NA, NSA = NA, TEP = "Y", AMN = "Y", AMS = "Y",   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "CAPRIMULGIFORMES",
          FAMILY = "TROCHILIDAE",
          Species_Latin_Name = "Phaethornis longirostris", # Long-billed Hermit, formerly Phaethornis (superciliosus) superciliosus
          GENUS = "Phaethornis",
          SPECIES = "longirostris",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "H", # Sensitivity
          STRAT = "U", # Foraging Strata
          CNTAB = "LT", # Center of Abundance    
          REL = "C", # Relative Abundance
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 2500, # Maximum Elevation (Changed from 1400)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA, # Microhabitats
          NHAB = 6, # Number of Habitats
          NZOO = 4, # Number of Zoogeographic Regions (Chris added F15, F11, F8)
          HAB1 = "F1", HAB2 = "F15", HAB3 = "F4", HAB4 = "F7", HAB5 = "F11", HAB6 = "F8", HAB7 = NA, #Habitat
          F1 = "Y", F2 = NA, F3 = NA, F4 = "Y", F5 = NA, F6 = NA, F7 = "Y", F8 = "Y",
          F9 = NA, F10 = NA, F11 = "Y", F12 = NA, F13 = NA, F14 = NA, F15 = "Y",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = "Y", CDH = NA, CHO = "Y",    
          EPC = "Y", STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = "Y", TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "CAPRIMULGIFORMES",
          FAMILY = "TROCHILIDAE",
          Species_Latin_Name = "Phaethornis striigularis", # Stripe-throated Hermit, formerly Phaethornis longuemareus
          GENUS = "Phaethornis",
          SPECIES = "striigularis",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "M", # Sensitivity
          STRAT = "U", # Foraging Strata
          CNTAB = "LT", # Center of Abundance    
          REL = "C", # Relative Abundance
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1570, # Maximum Elevation (Changed from 1200)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = "T", # Microhabitats 
          NHAB = 4, # Number of Habitats
          NZOO = 6, # Number of Zoogeographic Regions (Chris added F7, F14, MAH, EPC, removed AMN, AMS)
          HAB1 = "F1", HAB2 = "F15", HAB3 = "F7", HAB4 = "N14", HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "Y", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = "Y", F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "Y",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = "Y",
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = "Y", MAH = "Y", GCS = "Y", CDH = NA, CHO = "Y",    
          EPC = "Y", STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = "Y", TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "TITYRIDAE",
          Species_Latin_Name = "Schiffornis veraepacis", # Northern Schiffornis, formerly Schiffornis turdinus
          GENUS = "Schiffornis",
          SPECIES = "veraepacis",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "M", # Sensitivity
          STRAT = "U", # Foraging Strata
          CNTAB = "LT", # Center of Abundance (Chris changed to LT)
          REL = "U/F", # Relative Abundance (Chris added U)
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1700, # Maximum Elevation
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA, # Microhabitats 
          NHAB = 3, # Number of Habitats
          NZOO = 4, # Number of Zoogeographic Regions (Chris added CDH, CHO, EPC, removed ATL, added F15)
          HAB1 = "F1", HAB2 = "F4", HAB3 = "F15", HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "Y", F2 = NA, F3 = NA, F4 = "Y", F5 = NA, F6 = NA, F7 = NA, F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "Y",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = "Y", CDH = "Y", CHO = "Y",    
          EPC = "Y", STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "TITYRIDAE",
          Species_Latin_Name = "Schiffornis turdina", # Brown-winged Schiffornis, formerly Schiffornis turdinus
          GENUS = "Schiffornis",
          SPECIES = "turdina",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "H", # Sensitivity
          STRAT = "U", # Foraging Strata
          CNTAB = "LT", # Center of Abundance
          REL = "U/F", # Relative Abundance (Chris added U)
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1500, # Maximum Elevation
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 2, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA, # Microhabitats 
          NHAB = 3, # Number of Habitats
          NZOO = 4, # Number of Zoogeographic Regions (Chris removed GCS, CDH, CHO, added F15)
          HAB1 = "F1", HAB2 = "F4", HAB3 = "F15", HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "Y", F2 = NA, F3 = NA, F4 = "Y", F5 = NA, F6 = NA, F7 = NA, F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "Y",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = "Y", TEP = NA, AMN = "Y", AMS = "Y",   
          CSA = NA, ATL = "Y", PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "THRAUPIDAE",
          Species_Latin_Name = "Sporophila corvina", # Variable Seedeater, formerly Sporophila (americana) aurita
          GENUS = "Sporophila",
          SPECIES = "corvina",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "L", # Sensitivity
          STRAT = "U/M", # Foraging Strata
          CNTAB = "LT", # Center of Abundance
          REL = "F/C", # Relative Abundance (Chris added F)
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1500, # Maximum Elevation
          QMAX = "+", # Maximum elevation qualifier (Chris added +)
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA, # Microhabitats 
          NHAB = 5, # Number of Habitats
          NZOO = 3, # Number of Zoogeographic Regions (Chris added E to F7 and F15, and N13)
          HAB1 = "N14", HAB2 = "F1E", HAB3 = "F7E", HAB4 = "F15E", HAB5 = "N13", HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "E", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = "E", F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "E",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = "Y", N14 = "Y",
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = "Y", CDH = NA, CHO = "Y",    
          EPC = "Y", STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "THAMNOPHILIDAE",
          Species_Latin_Name = "Thamnophilus atrinucha", # Black-crowned Antshrike, formerly Thamnophilus punctatus
          GENUS = "Thamnophilus",
          SPECIES = "atrinucha",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "L", # Sensitivity
          STRAT = "U/M", # Foraging Strata
          CNTAB = "LT", # Center of Abundance
          REL = "C", # Relative Abundance
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1500, # Maximum Elevation
          QMAX = NA, # Maximum elevation qualifier (Chris added +)
          CP = 4, # Conservation Priority
          RP = 2, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA, # Microhabitats 
          NHAB = 3, # Number of Habitats
          NZOO = 3, # Number of Zoogeographic Regions (Chris removed F12 and F8)
          HAB1 = "F1E", HAB2 = "F7", HAB3 = "F15", HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "E", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = "Y", F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "Y",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = "Y", CDH = NA, CHO = "Y",    
          EPC = "Y", STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "THAMNOPHILIDAE",
          Species_Latin_Name = "Willisornis poecilinotus", # Common Scale-backed Antbird, formerly Hylophylax poecilinota
          GENUS = "Willisornis",
          SPECIES = "poecilinotus",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "M", # Sensitivity
          STRAT = "U", # Foraging Strata
          CNTAB = "LT", # Center of Abundance
          REL = "F", # Relative Abundance
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1350, # Maximum Elevation (Chris changed from 950)
          QMAX = "+", # Maximum elevation qualifier (Chris added +)
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA, # Microhabitats 
          NHAB = 1, # Number of Habitats
          NZOO = 2, # Number of Zoogeographic Regions
          HAB1 = "F1", HAB2 = NA, HAB3 = NA, HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "Y", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = NA, F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = NA,
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = "Y", AMS = "Y",   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "FURNARIIDAE",
          Species_Latin_Name = "Xiphorhynchus susurrans", # Cocoa Woodcreeper, formerly Xiphorhynchus (guttatus) guttatus
          GENUS = "Xiphorhynchus",
          SPECIES = "susurrans",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "L", # Sensitivity
          STRAT = "U/C", # Foraging Strata
          CNTAB = "LT", # Center of Abundance
          REL = "F/C", # Relative Abundance (Chris added F)
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 2400, # Maximum Elevation (Chris changed from 1100)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA, # Microhabitats 
          NHAB = 4, # Number of Habitats
          NZOO = 5, # Number of Zoogeographic Regions (Chris removed F2, AMN, AMS, ATL; added F8 and F7)
          HAB1 = "F1", HAB2 = "F8", HAB3 = "F7", HAB4 = "F15", HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "Y", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = "Y", F8 = "Y",
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = NA,
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = "Y", CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = "Y", TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "FURNARIIDAE",
          Species_Latin_Name = "Xiphorhynchus guttatus", # Buff-throated Woodcreeper, formerly Xiphorhynchus (guttatus) guttatus
          GENUS = "Xiphorhynchus",
          SPECIES = "guttatus",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "L", # Sensitivity
          STRAT = "U/C", # Foraging Strata
          CNTAB = "LT", # Center of Abundance
          REL = "F/C", # Relative Abundance (Chris added F)
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1250, # Maximum Elevation (Chris changed from 1100)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA, # Microhabitats 
          NHAB = 7, # Number of Habitats
          NZOO = 5, # Number of Zoogeographic Regions (Chris added F7, F8, F12, F13, removed GCS; added TEP, ATL)
          HAB1 = "F1", HAB2 = "F2", HAB3 = "F15", HAB4 = "F7", HAB5 = "F8", HAB6 = "F12", HAB7 = "F13", #Habitat
          F1 = "Y", F2 = "Y", F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = "Y", F8 = "Y",
          F9 = NA, F10 = NA, F11 = NA, F12 = "Y", F13 = "Y", F14 = NA, F15 = "Y",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = "Y", TEP = "Y", AMN = "Y", AMS = "Y",   
          CSA = NA, ATL = "Y", PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "CUCULIFORMES",
          FAMILY = "CUCULIDAE",
          Species_Latin_Name = "Coccyzus longirostris", # Hispaniolan Lizard-Cuckoo, Caribbean species never included in Parker
          GENUS = "Coccyzus",
          SPECIES = "longirostris",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "L", # Sensitivity (Chris is unsure of what to log this as)
          STRAT = "M/C", # Foraging Strata
          CNTAB = "LT", # Center of Abundance (Chris is unsure of what to log this as, they seem pretty widespread)
          REL = "F/C", # Relative Abundance 
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 2200, # Maximum Elevation
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority 
          RP = 2, # Research Priority (Chris is unsure of what to log this as)
          SUB = "HIS", # Subregions
          MICRO = NA, # Microhabitats 
          NHAB = 3, # Number of Habitats (Chris is unsure of N2)
          NZOO = 1, # Number of Zoogeographic Regions
          HAB1 = "F7", HAB2 = "F1", HAB3 = "N2", HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "Y", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = "Y", F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = NA,
          N1 = NA, N2 = "Y", N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = "Y", LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "TYRANNIDAE",
          Species_Latin_Name = "Contopus hispaniolensis", # Hispaniolan Pewee, Caribbean species never included in Parker
          GENUS = "Contopus",
          SPECIES = "hispaniolensis",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "L", # Sensitivity (Chris is unsure of what to log this as)
          STRAT = "U", # Foraging Strata
          CNTAB = "LT", # Center of Abundance (Chris is unsure of what to log this as, they seem pretty widespread)
          REL = "C", # Relative Abundance 
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1800, # Maximum Elevation
          QMAX = "+", # Maximum elevation qualifier
          CP = 4, # Conservation Priority 
          RP = 2, # Research Priority (Chris is unsure of what to log this as)
          SUB = "HIS", # Subregions
          MICRO = NA, # Microhabitats 
          NHAB = 4, # Number of Habitats (Chris is unsure of these, not much information on BOW)
          NZOO = 1, # Number of Zoogeographic Regions
          HAB1 = "F10E", HAB2 = "F7E", HAB3 = "F1E", HAB4 = "F15E", HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "E", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = "E", F8 = NA,
          F9 = NA, F10 = "E", F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "E",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = "Y", LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "THRAUPIDAE",
          Species_Latin_Name = "Melopyrrha portoricensis", # Puerto Rican Bullfinch, Caribbean species never included in Parker
          GENUS = "Melopyrrha",
          SPECIES = "portoricensis",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "L", # Sensitivity
          STRAT = "T/U/M/C", # Foraging Strata
          CNTAB = "MM", # Center of Abundance (Chris is unsure of what to log this as, they seem pretty widespread)
          REL = "C", # Relative Abundance 
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1000, # Maximum Elevation (This information was not present on BOW as of 12/2/20)
          QMAX = "+", # Maximum elevation qualifier
          CP = 4, # Conservation Priority 
          RP = 3, # Research Priority (Chris is unsure of what to log this as)
          SUB = "PUE", # Subregions
          MICRO = NA, # Microhabitats 
          NHAB = 3, # Number of Habitats (Chris is unsure of these, not much information on BOW)
          NZOO = 1, # Number of Zoogeographic Regions
          HAB1 = "F4", HAB2 = "N1", HAB3 = "F14", HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = NA, F2 = NA, F3 = NA, F4 = "Y", F5 = NA, F6 = NA, F7 = NA, F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = "Y", F15 = NA,
          N1 = "Y", N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = "Y", LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "SPINDALIDAE",
          Species_Latin_Name = "Spindalis portoricensis", # Puerto Rican Spindalis, Caribbean species never included in Parker
          GENUS = "Spindalis",
          SPECIES = "portoricensis",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "L", # Sensitivity
          STRAT = "T/U/M/C", # Foraging Strata
          CNTAB = "LT", # Center of Abundance (Chris is unsure of what to log this as, they seem pretty widespread)
          REL = "C", # Relative Abundance 
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1000, # Maximum Elevation (This information was not present on BOW as of 12/2/20)
          QMAX = "+", # Maximum elevation qualifier
          CP = 4, # Conservation Priority 
          RP = 3, # Research Priority (Chris is unsure of what to log this as)
          SUB = "PUE", # Subregions
          MICRO = NA, # Microhabitats 
          NHAB = 2, # Number of Habitats (Chris is unsure of these, not much information on BOW)
          NZOO = 1, # Number of Zoogeographic Regions
          HAB1 = "F15E", HAB2 = "F1E", HAB3 = NA, HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "E", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = NA, F8 = "Y",
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "E",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = "Y", LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = NA, AMS = NA,   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "THAMNOPHILIDAE",
          Species_Latin_Name = "Epinecrophylla haematonota", # Rufous-backed Stipplethroat, formerly known as Myrmotherula (haematonota) haematonota
          GENUS = "Epinecrophylla",
          SPECIES = "haematonota",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "H", # Sensitivity
          STRAT = "U", # Foraging Strata
          CNTAB = "LT", # Center of Abundance    
          REL = "U/F", # Relative Abundance
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1300, # Maximum Elevation (Chris changed to 1300)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 2, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA, # Microhabitats
          NHAB = 1, # Number of Habitats
          NZOO = 2, # Number of Zoogeographic Regions
          HAB1 = "F1", HAB2 = NA, HAB3 = NA, HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "Y", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = NA, F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = NA,
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = "Y", AMS = "Y",   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "TYRANNIDAE",
          Species_Latin_Name = "Hemitriccus griseipectus", # White-bellied Tody-Tyrant, formerly known as Hemitriccus zosterops
          GENUS = "Hemitriccus",
          SPECIES = "griseipectus",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "L", # Sensitivity (Chris changed to L)
          STRAT = "M", # Foraging Strata
          CNTAB = "LT", # Center of Abundance    
          REL = "C", # Relative Abundance (Chris changed to C)
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 900, # Maximum Elevation (Chris changed to 900)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 2, # Research Priority
          SUB = NA, # Subregions
          MICRO = "V", # Microhabitats (Chris added V)
          NHAB = 1, # Number of Habitats
          NZOO = 2, # Number of Zoogeographic Regions
          HAB1 = "F1", HAB2 = NA, HAB3 = NA, HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "Y", F2 = NA, F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = NA, F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = NA,
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = NA, AMS = "Y",   
          CSA = NA, ATL = "Y", PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "CAPITONIDAE",
          Species_Latin_Name = "Capito auratus", # Gilded Barbet, formerly known as Capito niger
          GENUS = "Capito",
          SPECIES = "auratus",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "M", # Sensitivity
          STRAT = "M/C", # Foraging Strata
          CNTAB = "LT", # Center of Abundance    
          REL = "C", # Relative Abundance
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 1700, # Maximum Elevation (Chris changed to 1700)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA,
          NHAB = 3, # Number of Habitats
          NZOO = 2, # Number of Zoogeographic Regions
          HAB1 = "F1", HAB2 = "F15", HAB3 = "F2", HAB4 = NA, HAB5 = NA, HAB6 = NA, HAB7 = NA, #Habitat
          F1 = "Y", F2 = "Y", F3 = NA, F4 = NA, F5 = NA, F6 = NA, F7 = NA, F8 = NA,
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = "Y",
          N1 = NA, N2 = NA, N3 = NA, N4 = NA, N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = NA,
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = "Y", AMS = "Y",   
          CSA = NA, ATL = NA, PAM = NA, PAT = NA,
          STATUS = NA) %>% # resident
  add_row(ORDER = "PASSERIFORMES",
          FAMILY = "TYRANNIDAE",
          Species_Latin_Name = "Nesotriccus murina", # Southern Mouse-colored Tyrannulet, formerly known as Phaeomyias murina
          GENUS = "Nesotriccus",
          SPECIES = "murina",
          NUMB = NA, # Parker classification number, irrelevant
          SNST = "L", # Sensitivity
          STRAT = "M", # Foraging Strata (Chris changed to M)
          CNTAB = "LT", # Center of Abundance    
          REL = "F/P", # Relative Abundance
          MIN = 0, # Minimum Elevation, 0 = Lowlands
          QMIN = NA, # Minimum elevation qualifier, ? = Uncertain value.
          MAX = 2400, # Maximum Elevation (Chris changed to 2400)
          QMAX = NA, # Maximum elevation qualifier
          CP = 4, # Conservation Priority
          RP = 3, # Research Priority
          SUB = NA, # Subregions
          MICRO = NA,
          NHAB = 8, # Number of Habitats (Chris added N2, arid montane scrub, and removed F14, mangrove forest)
          NZOO = 4, # Number of Zoogeographic Regions
          HAB1 = "N1", HAB2 = "N2", HAB3 = "N4", HAB4 = "N14", HAB5 = "F7", HAB6 = "F8", HAB7 = "F3", #Habitat
          F1 = NA, F2 = NA, F3 = "Y", F4 = NA, F5 = NA, F6 = NA, F7 = "Y", F8 = "Y",
          F9 = NA, F10 = NA, F11 = NA, F12 = NA, F13 = NA, F14 = NA, F15 = NA,
          N1 = "Y", N2 = "Y", N3 = NA, N4 = "Y", N5 = NA, N6 = NA, N7 = NA, N8 = NA,
          N9 = NA, N10 = NA, N11 = NA, N12 = NA, N13 = NA, N14 = "Y",
          A1 = NA, A2 = NA, A3 = NA, A4 = NA, A5 = NA, A6 = NA, A7 = NA, A8 = NA,
          A9 = NA, A10 = NA, A11 = NA, A12 = NA,
          GAN = NA, LAN = NA, BSR = NA, MPL = NA, PAS = NA, MAH = NA, GCS = NA, CDH = NA, CHO = NA,    
          EPC = NA, STP = NA, NAN = NA, CAN = NA, SAN = NA, NSA = NA, TEP = NA, AMN = "Y", AMS = "Y",   
          CSA = "Y", ATL = "Y", PAM = NA, PAT = NA,
          STATUS = NA) # resident

# Which species were STILL not joined due to taxonomic changes or subspecies? — eventually will be 0
unjoined <- anti_join(TRACEData, parker, by = "Species_Latin_Name")
unique(unjoined$Species_Latin_Name)

# Missing species that still need to be added:
# This not an immediate priority as of 1/10/22 because these species don't have Hg samples
# Lonchura malacca - Tricolored Munia
# Psittacara chloropterus - Hispaniolan Parakeet

# Taxa that cannot be joined (this is okay):
# Aves sp.
# Hemitriccus sp.- Tody-Tyrants
# Ramphotrigon sp. - Flatbills
# Xiphorhynchus sp. - Woodcreepers
# Myrmotherula sp. - Antwrens


##### CLASSIFYING SPECIES BY TROPHIC NICHE USING PIGOT et al. (2020) CRITERIA #####

# Reading in Supplementary Dataset 1 from Pigot et al. (2020), which features trophic and foraging niches
pigot <- read_excel("Pigot2020TrophicNiche.xlsx") %>%
  clean_names(case = "mixed") %>% # cleaning up column names
  mutate(Species_Latin_Name = str_replace(Binomial, "_", " ")) %>% # getting rid of underscores so we can merge datasets better
  select(Species_Latin_Name, Trophic_Level, Trophic_Niche, Foraging_Niche)

# Which species were not joined due to taxonomic changes or subspecies since 2019?
unjoined <- anti_join(TRACEData, pigot, by = "Species_Latin_Name")
unique(unjoined$Species_Latin_Name) # We have 67 taxa that were unjoined

# New name is being assigned to old Pigot name
pigot <- pigot %>%
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Thryothorus modestus", "Cantorchilus modestus", Species_Latin_Name)) %>% # Cabanis's Wren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Thryothorus semibadius", "Cantorchilus semibadius", Species_Latin_Name)) %>% # Riverside Wren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Thryothorus thoracicus", "Cantorchilus thoracicus", Species_Latin_Name)) %>% # Stripe-breasted Wren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Wilsonia pusilla", "Cardellina pusilla", Species_Latin_Name)) %>% # Wilson's Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Pipra mentalis", "Ceratopipra mentalis", Species_Latin_Name)) %>% # Red-capped Manakin
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Cercomacra tyrannina", "Cercomacroides tyrannina", Species_Latin_Name)) %>% # Dusky Antbird
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Chlorospingus ophthalmicus", "Chlorospingus flavopectus", Species_Latin_Name)) %>% # Common Chlorospingus
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Oporornis formosus", "Geothlypis formosa", Species_Latin_Name)) %>% # Kentucky Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Oporornis philadelphia", "Geothlypis philadelphia", Species_Latin_Name)) %>% # Mourning Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Oporornis philadelphia", "Geothlypis philadelphia", Species_Latin_Name)) %>% # Mourning Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Oporornis tolmiei", "Geothlypis tolmiei", Species_Latin_Name)) %>% # MacGillivray's Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Myrmeciza fortis", "Hafferia fortis", Species_Latin_Name)) %>% # Sooty Antbird
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Myrmotherula hauxwelli", "Isleria hauxwelli", Species_Latin_Name)) %>% # Plain-throated Antwren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Vermivora peregrina", "Leiothlypis peregrina", Species_Latin_Name)) %>% # Tennessee Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Larus atricilla", "Leucophaeus atricilla", Species_Latin_Name)) %>% # Laughing Gull
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Loxigilla portoricensis", "Melopyrrha portoricensis", Species_Latin_Name)) %>% # Puerto Rican Bullfinch
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Schistocichla leucostigma", "Myrmelastes leucostigma", Species_Latin_Name)) %>% # Spot-winged Antbird
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Gymnopithys salvini", "Oneillornis salvini", Species_Latin_Name)) %>% # White-throated Antbird
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Parula gutturalis", "Oreothlypis gutturalis", Species_Latin_Name)) %>% # Flame-throated Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Seiurus motacilla", "Parkesia motacilla", Species_Latin_Name)) %>% # Lousiana Waterthrush
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Seiurus noveboracensis", "Parkesia noveboracensis", Species_Latin_Name)) %>% # Northern Waterthrush
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Thryothorus maculipectus", "Pheugopedius maculipectus", Species_Latin_Name)) %>% # Spot-breasted Wren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Cyanocorax morio",  "Psilorhinus morio", Species_Latin_Name)) %>% # Brown Jay
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Aratinga chloroptera", "Psittacara chloropterus", Species_Latin_Name)) %>% # Hispaniolan Parakeet
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Buteo magnirostris", "Rupornis magnirostris", Species_Latin_Name)) %>% # Roadside Hawk
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Parula americana", "Setophaga americana", Species_Latin_Name)) %>% # Northern Parula
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Wilsonia citrina", "Setophaga citrina", Species_Latin_Name)) %>% # Hooded Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Dendroica discolor", "Setophaga discolor", Species_Latin_Name)) %>% # Prairie Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Dendroica magnolia", "Setophaga magnolia", Species_Latin_Name)) %>% # Magnolia Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Dendroica pensylvanica", "Setophaga pensylvanica", Species_Latin_Name)) %>% # Chestnut-sided Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Dendroica petechia", "Setophaga petechia", Species_Latin_Name)) %>% # Yellow Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Dendroica tigrina", "Setophaga tigrina", Species_Latin_Name)) %>% # Cape May Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Oryzoborus angolensis", "Sporophila angolensis", Species_Latin_Name)) %>% # Chestnut-bellied Seed-Finch
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Oryzoborus funereus", "Sporophila funerea", Species_Latin_Name)) %>% # Thick-billed Seed-Finch
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Sporophila torqueola", "Sporophila morelleti", Species_Latin_Name)) %>% # Morelet's Seedeater
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Tangara larvata", "Stilpnia larvata", Species_Latin_Name)) %>% # Golden-hooded Tanager
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Thryothorus rufalbus", "Thryophilus rufalbus", Species_Latin_Name)) %>% # Rufous-and-white Wren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Vermivora pinus", "Vermivora cyanoptera", Species_Latin_Name)) %>% # Blue-winged Warbler
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Myrmeciza hemimelaena", "Sciaphylax hemimelaena", Species_Latin_Name)) %>% # Chestnut-tailed Antbird (this  species was split but then lumped back)
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Myrmeciza hyperythra", "Myrmelastes hyperythrus", Species_Latin_Name)) %>% # Plumbeous Antbird
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Myrmeciza melanoceps", "Akletos melanoceps", Species_Latin_Name)) %>%  # White-shouldered Antbird
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Myrmeciza atrothorax", "Myrmophylax atrothorax", Species_Latin_Name)) %>% # Black-throated Antbird
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Veniliornis passerinus", "Dryobates passerinus", Species_Latin_Name)) %>% # Little Woodpecker
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Thryothorus genibarbis", "Pheugopedius genibarbis", Species_Latin_Name)) %>% # Moustached Wren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Oryzoborus atrirostris", "Sporophila atrirostris", Species_Latin_Name)) %>% # Black-billed Seed-Finch
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Pipra chloromeros", "Ceratopipra chloromeros", Species_Latin_Name)) %>% # Round-tailed Manakin
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Thryothorus leucotis", "Cantorchilus leucotis", Species_Latin_Name)) %>% # Buff-breasted Wren
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Amazilia candida", "Chlorestes candida", Species_Latin_Name)) %>% # White-bellied Emerald
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Chlorostilbon canivetii", "Cynanthus canivetii", Species_Latin_Name)) %>% # Canivet's Emerald
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Glaucidium brasilianum", "Nannopterum brasilianum", Species_Latin_Name)) %>% # Neotropical Cormorant
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Hylocharis cyanus", "Chlorestes cyanus", Species_Latin_Name)) %>% # White-chinned Sapphire
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Pteroglossus beauharnaesii", "Pteroglossus beauharnaisii", Species_Latin_Name)) %>% # Curl-crested Aracari
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Tachyphonus luctuosus", "Loriotus luctuosus", Species_Latin_Name)) %>% # White-shouldered Tanager
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Amazilia lactea", "Chionomesa lactea", Species_Latin_Name)) %>% # Sapphire-spangled Emerald
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Picoides fumigatus", "Dryobates fumigatus", Species_Latin_Name)) %>% # Smoky-brown Woodpecker
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Icterus (dominicensis) prosthemelas", "Icterus prosthemelas", Species_Latin_Name)) %>% # Black-cowled Oriole
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Hylophilus decurtatus", "Pachysylvia decurtata", Species_Latin_Name)) %>% # Lesser Greenlet
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Leucopternis schistaceus", "Buteogallus schistaceus", Species_Latin_Name)) %>% # Slate-colored Hawk
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Speotyto cunicularia", "Athene cunicularia", Species_Latin_Name)) %>% # Burrowing Owl
  mutate(Species_Latin_Name = ifelse(Species_Latin_Name == "Phaeomyias murina", "Nesotriccus murina", Species_Latin_Name)) # Southern Mouse-colored Tyrannulet

# Which species were STILL not joined due to taxonomic changes or subspecies? — eventually will be 0
unjoined <- anti_join(TRACEData, pigot, by = "Species_Latin_Name")
unique(unjoined$Species_Latin_Name) # We have 10 taxa that have been split/lumped since 2019 or were just excluded from the df

# Manually adding species to accommodate for splits/lumps information
# I am choosing to leave the original classifications in the df so that we can reference it later if necessary
# For splits, each species gets a new row with the same trophic information defined by Pigot et al. (2020)

pigot <- pigot %>% 
  add_row(Species_Latin_Name = "Cyanoloxia cyanoides",
          Trophic_Level = "Herbivore",
          Trophic_Niche = "Granivore", # Blue-Black Grosbeak, formerly Cyanocompsa cyanoides
          Foraging_Niche = "Granivore arboreal") %>%
  add_row(Species_Latin_Name = "Cyanoloxia rothschildii",
          Trophic_Level = "Herbivore",
          Trophic_Niche = "Granivore", # Amazonian Grosbeak, formerly Cyanocompsa cyanoides
          Foraging_Niche = "Granivore arboreal") %>% 
  add_row(Species_Latin_Name = "Schiffornis veraepacis",
          Trophic_Level = "Omnivore",
          Trophic_Niche = "Omnivore", # Northern Schiffornis, formerly Schiffornis turdina
          Foraging_Niche = "NA") %>% 
  add_row(Species_Latin_Name = "Momotus lessonii",
          Trophic_Level = "Carnivore",
          Trophic_Niche = "Omnivore", # Lesson's Motmot, formerly Momotus momota
          Foraging_Niche = "Generalist") %>%
  add_row(Species_Latin_Name = "Epinecrophylla amazonica",
          Trophic_Level = "Carnivore",
          Trophic_Niche = "Invertivore", # Rio Madeira Stipplethroat, formerly Epinecrophylla haematonota
          Foraging_Niche = "Invertivore glean arboreal") %>% 
  add_row(Species_Latin_Name = "Hemitriccus sp.",
          Trophic_Level = "Carnivore",
          Trophic_Niche = "Invertivore", # Tody-Tyrant sp.
          Foraging_Niche = "Invertivore sally surface") %>% 
  add_row(Species_Latin_Name = "Ramphotrigon sp.",
          Trophic_Level = "Carnivore",
          Trophic_Niche = "Invertivore", # Flatbill sp.
          Foraging_Niche = "Invertivore sally surface") %>%
  add_row(Species_Latin_Name = "Myrmotherula sp.",
          Trophic_Level = "Carnivore",
          Trophic_Niche = "Invertivore", # Antwren sp.
          Foraging_Niche = "Invertivore glean arboreal") %>% 
  add_row(Species_Latin_Name = "Xiphorhynchus sp.",
          Trophic_Level = "Carnivore",
          Trophic_Niche = "Invertivore", # Woodcreeper sp.
          Foraging_Niche = "Invertivore bark") %>% 
  mutate(Foraging_Niche = ifelse(Foraging_Niche == "NA", "Generalist", Foraging_Niche))
# NOTE: Omnivores do not have foraging niches in this df, but I am doing this for visualization purposes later

# For added ecological clarity, we change trophic niche "Vertivore" to "Terrestrial vertivore"
pigot <- pigot %>% 
  mutate(Trophic_Niche = ifelse(Trophic_Niche == "Vertivore", "Terrestrial vertivore", Trophic_Niche))

# Which species were STILL not joined due to taxonomic changes or subspecies? — eventually will be 0
unjoined <- anti_join(TRACEData, pigot, by = "Species_Latin_Name")
unique(unjoined$Species_Latin_Name) # Should be 0 — it is!


# JOINING DATASETS --------------------------------------------

# Phew! That was a lot of work. Let's put everything together now.
CollectiveData <- left_join(TRACEData, taxa, by = "Species_Latin_Name") %>%
  left_join(pigot, by = "Species_Latin_Name") %>%
  left_join(parker, by = "Species_Latin_Name") %>%
  unite("Family", c(Family.x, Family.y), na.rm = T, remove = T, sep = "") %>%
  # mending the unite() weirdness
  mutate(Family = na_if(Family, "")) %>% 
  select(-c(ORDER:NUMB)) %>% # Getting rid of superfluous information
  # simplifying primary habitat association for modeling purposes
  # Aquatic, Forest, and Grassland/scrub
  #mutate(Primary_Habitat = ifelse(HAB1 %in% c("F1", "F1E", "F2", "F3", "F3?", "F4", "F4E", "F5", "F6", "F7", "F7E", "F8",
  #                                            "F9", "F10", "F10E", "F11", "F12", "F13", "F14", "F15", "F15E"), "Forest",
  #                                                     ifelse(HAB1 %in% c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8",
  #                                                                        "N9", "N10", "N11", "N12", "N13", "N14"), "Grassland/scrub",
  #                                                            ifelse(HAB1 %in% c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8",
  #                                                                               "A9", "A10", "A11", "A12"), "Aquatic",
  #                                                                   HAB1)))) %>% 
  mutate(Primary_Habitat = ifelse(HAB1 %in% c("F1", "F1E", "F2", "F10",
                                             "F10E", "F11", "F12", "F13", "F14"), "Lowland evergreen forest",
                                  ifelse(HAB1 %in% c("F4", "F4E", "F5", "F6"), "Montane evergreen forest",
                                         ifelse(HAB1 %in% c("F7", "F7E", "F9"), "Lowland deciduous forest",
                                                ifelse(HAB1 %in% c("F3", "F3?", "F8", "F15", "F15E"), "Secondary forest", # river-edge, gallery, and true secondary
                                                       ifelse(HAB1 %in% c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8",
                                                                          "N9", "N10", "N11", "N12", "N13", "N14"), "Grassland/scrub",
                                                              ifelse(HAB1 %in% c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8",
                                                                                 "A9", "A10", "A11", "A12"), "Aquatic",
                                                                     HAB1))))))) %>%
  # Now let's make the species habitat types more readable
  mutate(HAB1 = ifelse(HAB1 %in% c("F1", "F1E"), "Tropical lowland evergreen forest", HAB1)) %>% 
  mutate(HAB1 = ifelse(HAB1 == "F2", "Flooded tropical evergreen forest", HAB1)) %>% 
  mutate(HAB1 = ifelse(HAB1 %in% c("F3", "F3?"), "River-edge forest", HAB1)) %>% 
  mutate(HAB1 = ifelse(HAB1 %in% c("F4", "F4E"), "Montane evergreen forest", HAB1)) %>% 
  mutate(HAB1 = ifelse(HAB1 == "F5", "Elfin forest", HAB1)) %>% 
  mutate(HAB1 = ifelse(HAB1 == "F6", "Polylepis woodland", HAB1)) %>% 
  mutate(HAB1 = ifelse(HAB1 %in% c("F7", "F7E"), "Tropical deciduous forest", HAB1)) %>% 
  mutate(HAB1 = ifelse(HAB1 == "F8", "Gallery forest", HAB1)) %>% 
  mutate(HAB1 = ifelse(HAB1 == "F9", "Southern temperate forest", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 %in% c("F10", "F10E"), "Pine forest", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "F11", "Pine-Oak forest", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "F12", "White sand forest", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "F13", "Palm forest", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "F14", "Mangrove forest", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 %in% c("F15", "F15E"), "Secondary forest", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "N1", "Arid lowland scrub", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "N2", "Arid montane scrub", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "N3", "Semihumid/humid montane scrub", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "N4", "Cerrado", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "N5", "Campo grasslands", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "N6", "Low, seasonally wet grassland", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "N7", "Southern temperate grassland", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "N8", "Northern temperate grassland", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "N9", "Puna grassland", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "N10", "Paramo grassland", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "N11", "Riparian thickets", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "N12", "River island scrub", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "N13", "Pastures/agricultural lands", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "N14", "Second-growth scrub", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "A1", "Freshwater marshes", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "A2", "Saltwater/brackish marshes", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "A3", "Coastal sand beaches/mudflats", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "A4", "Coastal rocky beaches", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "A5", "Riverine sand beaches", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "A6", "Freshwater lakes and ponds", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "A7", "Alkaline lakes", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "A8", "Rivers", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "A9", "Streams", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "A10", "Bogs", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "A11", "Coastal waters", HAB1)) %>%
  mutate(HAB1 = ifelse(HAB1 == "A12", "Pelagic waters", HAB1)) %>%
  # setting all the undefined species (e.g. Xiphorhynchus sp.) to "resident"
  mutate(Migratory_Status = if_else(is.na(Migratory_Status), "Resident", Migratory_Status)) %>% 
  # excluding feather samples from migratory species because of molt uncertainty
  mutate(Tail_Hg_ppm = ifelse(Migratory_Status != "Resident" & !is.na(Tail_Hg_ppm), NA, Tail_Hg_ppm)) %>%
  mutate(Body_Hg_ppm = ifelse(Migratory_Status != "Resident" & !is.na(Body_Hg_ppm), NA, Body_Hg_ppm)) %>%
  # renaming sites per CINCIA's temporary request
  mutate(Site_Name = if_else(Site_Name == "Azul Mine", "Mining Site A", Site_Name)) %>%
  mutate(Site_Name = if_else(Site_Name == "La Torre", "Reserva Nacional Tambopata", Site_Name)) %>%
  mutate(Site_Name = if_else(Site_Name == "Paolita Mine", "Mining Site B", Site_Name)) %>%
  mutate(Site_Name = if_else(Site_Name == "Santa Rita Mine", "Mining Site C", Site_Name))


# PRODUCING SUMMARY STATISTICS --------------------------------------------

# How many samples do we have and from what organizations
# this total includes multiple tissue samples from the same individual
CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  select(Collecting_Organization, Tissue_Type, Concentration) %>% 
  filter(!is.na(Concentration)) %>%
  count(Collecting_Organization) %>% 
  ungroup() %>% 
  mutate(Percent = (n/sum(n))*100) %>% 
  view()

# How many total sampled orders do we have?
summary <- CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  select(Order, Tissue_Type, Concentration) %>% 
  filter(!is.na(Concentration)) %>%
  count(Order) %>% 
  view()

# How many total samples do we have?
# this total includes multiple tissue samples from the same individual
sum(summary$n)

# How many total sampled families do we have?
summary <- CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  select(Family, Tissue_Type, Concentration) %>% 
  filter(!is.na(Concentration), !is.na(Family)) %>%
  count(Family) %>% 
  view()
nrow(summary)

# How many total sampled taxa do we have?
summary <- CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  select(Species_Common_Name, Tissue_Type, Concentration) %>% 
  filter(!is.na(Concentration)) %>%
  # only including full species names
  filter(!(str_detect(Species_Common_Name, " sp."))) %>% 
  count(Species_Common_Name) %>% 
  view()
nrow(summary)

# How many total sampled sites do we have?
summary <- CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  select(Site_Name, Tissue_Type, Concentration) %>% 
  filter(!is.na(Concentration), !is.na(Site_Name)) %>%
  count(Site_Name) %>% 
  view()
nrow(summary)

# How many total sampled stations do we have?
summary <- CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  select(Banding_Station_Name, Tissue_Type, Concentration) %>% 
  filter(!is.na(Concentration), !is.na(Banding_Station_Name)) %>%
  count(Banding_Station_Name) %>% 
  view()
nrow(summary)

# How many samples do we have and for what tissues?
summary <- CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  select(Country, Year, Order, Family, Species_Common_Name, Tissue_Type, Concentration) %>% 
  filter(!is.na(Concentration)) %>% 
  group_by(Tissue_Type) %>% 
  summarize(Country = unique(Country), Min_Year = min(Year), Max_Year = max(Year),
            n_Species = length(unique(Species_Common_Name)), n = n(),
            Mean = mean(Concentration), SD = sd(Concentration), Min = min(Concentration),
            Max = max(Concentration), CV = (SD/Mean)*100) %>% 
  mutate_at(5:11, funs(round(., 3))) %>%
  mutate(Year_Range = if_else(Min_Year != Max_Year,
                              str_c(Min_Year, "—", Max_Year), as.character(Max_Year)),
         Mean = if_else(n > 1, str_c(Mean, " ± ", SD), as.character(Mean)),
         Range = if_else(n > 1, str_c(Min, "—", Max), "NA")) %>% 
  select(-c(Min_Year, Max_Year, SD, Min, Max)) %>% 
  ungroup() %>% 
  mutate(Percent = (n/2176)*100) %>%
  view()

# How many samples do we have and for what countries?
CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Body_Hg_ppm, Tail_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>%
  select(Country, Tissue_Type, Concentration) %>% 
  filter(!is.na(Concentration)) %>%
  count(Country) %>% 
  ungroup() %>% 
  mutate(Percent = (n/sum(n))*100) %>% 
  view()

# Creating a country Hg table displaying arithmetic mean and SD
CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Body_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  select(Tissue_Type, Country, Year, Species_Common_Name, Concentration) %>% 
  filter(!is.na(Concentration)) %>%
  group_by(Country, Tissue_Type) %>%
  summarize(Country = unique(Country), Min_Year = min(Year), Max_Year = max(Year),
            n_Species = length(unique(Species_Common_Name)), n = n(),
            Mean = mean(Concentration), SD = sd(Concentration), Min = min(Concentration),
            Max = max(Concentration), CV = (SD/Mean)*100) %>% 
  mutate_at(5:11, funs(round(., 3))) %>%
  mutate(Year_Range = if_else(Min_Year != Max_Year,
                              str_c(Min_Year, "—", Max_Year), as.character(Max_Year)),
         Mean = if_else(n > 1, str_c(Mean, " ± ", SD), as.character(Mean)),
         Range = if_else(n > 1, str_c(Min, "—", Max), "NA")) %>% 
  select(-c(Min_Year, Max_Year, SD, Min, Max)) %>% 
  view()

# Creating a site Hg table displaying arithmetic mean and SD
CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Body_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  select(Tissue_Type, Country, Site_Name, Year, Species_Common_Name, Concentration) %>% 
  filter(!is.na(Concentration)) %>%
  group_by(Country, Site_Name, Tissue_Type) %>%
  summarize(Site_Name = unique(Site_Name), Min_Year = min(Year), Max_Year = max(Year),
            n_Species = length(unique(Species_Common_Name)), n = n(),
            Mean = mean(Concentration), SD = sd(Concentration), Min = min(Concentration),
            Max = max(Concentration), CV = (SD/Mean)*100) %>% 
  mutate_at(5:12, funs(round(., 3))) %>%
  mutate(Year_Range = if_else(Min_Year != Max_Year,
                              str_c(Min_Year, "—", Max_Year), as.character(Max_Year)),
         Mean = if_else(n > 1, str_c(Mean, " ± ", SD), as.character(Mean)),
         Range = if_else(n > 1, str_c(Min, "—", Max), "NA")) %>% 
  select(-c(Min_Year, Max_Year, SD, Min, Max)) %>% 
  view()

# Creating an order Hg table displaying arithmetic mean and SD
CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Body_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  select(Order, Tissue_Type, Country, Year, Concentration) %>% 
  filter(!is.na(Concentration)) %>%
  group_by(Order, Tissue_Type) %>%
  summarize(Country = unique(Country), Min_Year = min(Year), Max_Year = max(Year), n = n(),
            Mean = mean(Concentration), SD = sd(Concentration), Min = min(Concentration),
            Max = max(Concentration), CV = (SD/Mean)*100) %>%
  mutate_at(7:11, funs(round(., 3))) %>%
  mutate(Year_Range = if_else(Min_Year != Max_Year,
                              str_c(Min_Year, "—", Max_Year), as.character(Max_Year)),
         Mean = if_else(n > 1, str_c(Mean, " ± ", SD), as.character(Mean)),
         Range = if_else(n > 1, str_c(Min, "—", Max), "NA")) %>% 
  select(-c(Min_Year, Max_Year, SD, Min, Max)) %>% 
  view() %>% 
  write.csv("Outputs/OrderHg.csv")

# Creating a family Hg table displaying arithmetic mean and SD
CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Body_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  select(Order, Family, Tissue_Type, Country, Year, Concentration) %>% 
  filter(!is.na(Concentration)) %>%
  group_by(Order, Family, Tissue_Type) %>%
  summarize(Country = unique(Country), Min_Year = min(Year), Max_Year = max(Year), n = n(),
            Mean = mean(Concentration), SD = sd(Concentration), Min = min(Concentration),
            Max = max(Concentration), CV = (SD/Mean)*100) %>% 
  mutate_at(7:12, funs(round(., 3))) %>%
  mutate(Year_Range = if_else(Min_Year != Max_Year,
                              str_c(Min_Year, "—", Max_Year), as.character(Max_Year)),
         Mean = if_else(n > 1, str_c(Mean, " ± ", SD), as.character(Mean)),
         Range = if_else(n > 1, str_c(Min, "—", Max), "NA")) %>% 
  select(-c(Min_Year, Max_Year, SD, Min, Max)) %>% 
  view() %>% 
  write.csv("Outputs/FamilyHg.csv")

# Creating a species Hg table displaying arithmetic mean and SD
CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Body_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  select(Order, Family, Species_Common_Name, Species_Latin_Name, Tissue_Type,
         Country, Year, Concentration) %>% 
  filter(!is.na(Concentration)) %>%
  group_by(Order, Family, Species_Common_Name, Species_Latin_Name, Tissue_Type) %>%
  summarize(Country = unique(Country), Min_Year = min(Year), Max_Year = max(Year), n = n(),
            Mean = mean(Concentration), SD = sd(Concentration), Min = min(Concentration),
            Max = max(Concentration), CV = (SD/Mean)*100) %>%
  mutate_at(9:14, funs(round(., 3))) %>% 
  mutate(Year_Range = if_else(Min_Year != Max_Year,
                              str_c(Min_Year, "—", Max_Year), as.character(Max_Year)),
         Mean = if_else(n > 1, str_c(Mean, " ± ", SD), as.character(Mean)),
         Range = if_else(n > 1, str_c(Min, "—", Max), "NA")) %>% 
  select(-c(Min_Year, Max_Year, SD, Min, Max)) %>% 
  view() %>% 
  write.csv("Outputs/SpeciesHg.csv")

# Creating a trophic niche Hg table displaying arithmetic mean and SD
CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Body_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  select(Trophic_Niche, Species_Common_Name, Tissue_Type, Concentration, Country, Year) %>% 
  filter(!is.na(Concentration)) %>%
  group_by(Trophic_Niche, Tissue_Type) %>%
  summarize(n_Species = length(unique(Species_Common_Name)), n = n(), Mean = mean(Concentration),
            SD = sd(Concentration), Min = min(Concentration), Max = max(Concentration), CV = (SD/Mean)*100,
            Country = unique(Country), Min_Year = min(Year), Max_Year = max(Year)) %>% 
  mutate_at(5:9, funs(round(., 3))) %>%
  mutate(Year_Range = if_else(Min_Year != Max_Year,
                              str_c(Min_Year, "—", Max_Year), as.character(Max_Year)),
         Mean = if_else(n > 1, str_c(Mean, " ± ", SD), as.character(Mean)),
         Range = if_else(n > 1, str_c(Min, "—", Max), "NA")) %>% 
  select(-c(Min_Year, Max_Year, SD, Min, Max)) %>% 
  view() %>% 
  write.csv("Outputs/TrophicHg.csv")

# Creating a primary habitat Hg table displaying arithmetic mean and SD
CollectiveData %>%
  pivot_longer(c(Blood_Hg_ppm, Tail_Hg_ppm, Body_Hg_ppm),
               names_to = "Tissue_Type", values_to = "Concentration") %>% 
  select(Primary_Habitat, Species_Common_Name, Tissue_Type, Concentration, Country, Year) %>% 
  filter(!is.na(Primary_Habitat), !is.na(Concentration)) %>%
  group_by(Primary_Habitat, Tissue_Type) %>%
  summarize(n_Species = length(unique(Species_Common_Name)), n = n(), Mean = mean(Concentration),
            SD = sd(Concentration), Min = min(Concentration), Max = max(Concentration), CV = (SD/Mean)*100,
            Country = unique(Country), Min_Year = min(Year), Max_Year = max(Year)) %>% 
  mutate_at(5:9, funs(round(., 3))) %>% 
  mutate(Year_Range = if_else(Min_Year != Max_Year,
                              str_c(Min_Year, "—", Max_Year), as.character(Max_Year)),
         Mean = if_else(n > 1, str_c(Mean, " ± ", SD), as.character(Mean)),
         Range = if_else(n > 1, str_c(Min, "—", Max), "NA")) %>% 
  select(-c(Min_Year, Max_Year, SD, Min, Max)) %>% 
  view() %>% 
  write.csv("Outputs/HabitatHg.csv")

unique(CollectiveData$HAB1)

# Creating species association table
CollectiveData %>%
  select(Order, Family, Species_Common_Name, Species_Latin_Name, Trophic_Level, Trophic_Niche,
         HAB1, Migratory_Status, Blood_Hg_ppm, Tail_Hg_ppm, Body_Hg_ppm) %>%
  # only including full species names
  filter(!(str_detect(Species_Common_Name, " sp."))) %>%
  filter(!is.na(Blood_Hg_ppm) | !is.na(Tail_Hg_ppm) | !is.na(Body_Hg_ppm)) %>%
  group_by(Order, Family, Species_Common_Name, Species_Latin_Name, Trophic_Level, Trophic_Niche,
           HAB1, Migratory_Status) %>%
  summarize(n = n()) %>% # Consolidating to 1 row per species
  view() %>% 
  write.csv("Outputs/SpeciesClass.csv")
