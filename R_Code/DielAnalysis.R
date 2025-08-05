# Diel activity

# Load some useful packages ----------------------------------------------------

library(lubridate)
library(ggplot2)
library(suncalc)
library(geosphere)
library(dplyr)
library(tidyr)
library(ggstatsplot)
library(lunar)

# set your working directory to your file folder -------------------------------

setwd("C:/Users/amajerus/Desktop/Annabelle/CSVs")

# Read in the data -------------------------------------------------------------

detections <- read.csv("detections_with_layers.csv") # this is a large file and will take some time to read in.
correction <- read.csv("camera_ID_corrected_for_duplicates.csv") # some locations have multiple IDs, this corrects it.
correction$X <- NULL
detections <- merge(detections, correction)

# create some useful columns ---------------------------------------------------

detections$DATE <- ymd(substr(detections$DETECTION_DATETIME, 1, 10))
detections$TIME <- substr(detections$DETECTION_DATETIME, 12, 19)
detections$TIME <- as.POSIXct(detections$TIME, format="%H:%M:%S")
detections$doy <- yday(detections$DATE)
OlsonNames() # these are all the valid time zone (tz) names for "getSunlightTimes"
detections$sunrise <- getSunlightTimes(data = data.frame(date=detections$DATE, lat=detections$LAT, lon=detections$LON), keep="sunrise", tz="Etc/GMT+5")$sunrise
detections$sunrise <- substr(detections$sunrise, 12, 19) # keep only time, not date
detections$sunrise <- as.POSIXct(detections$sunrise, format="%H:%M:%S") #make it a posix object
detections$sunset <- getSunlightTimes(data = data.frame(date=detections$DATE, lat=detections$LAT, lon=detections$LON), keep="sunset", tz="Etc/GMT+5")$sunset
detections$sunset <- substr(detections$sunset, 12, 19)
detections$sunset <- as.POSIXct(detections$sunset, format="%H:%M:%S")
detections$daylength <- daylength(lat=detections$LAT, doy=detections$doy) # hours of sunlight in the day
detections$dielCategory <- ifelse(detections$TIME>detections$sunrise & detections$TIME<detections$sunset, "DAY", "NIGHT") # did this detection happen while the sun was up?
detections$lunarIllumination <- lunar.illumination(detections$DATE)

#showing active time for species

SpeciesActiveModel = function(species) {
  species_data = subset(detections, detections$SPECIES==species)
  ggplot(species_data, aes(x=species_data$doy)) +
    geom_point(alpha=0.02, aes(y=as.POSIXct(species_data$TIME, format="%H:%M:%S"))) +
    geom_line(aes(y=as.POSIXct(species_data$sunrise, format="%H:%M:%S"), color="Sunrise")) +
    geom_line(aes(y=as.POSIXct(species_data$sunset, format="%H:%M:%S"), color="Sunset")) +
    scale_y_datetime(date_labels = "%H:%M") +
    ylab("Time") +
    xlab("Day of the year") +
    ggtitle(paste(species, "detections"))
  ggsave(paste0(species, "_Detections.png"), path = 'C:/Users/amajerus/Desktop/Annabelle/Figures/Detections/')
}

SpeciesActiveModel("Deer")
SpeciesActiveModel("Raccoon")
SpeciesActiveModel("Cottontail")
SpeciesActiveModel("Coyote")
SpeciesActiveModel("Other Bird")
SpeciesActiveModel("Squirrels and Chipmunks")
SpeciesActiveModel("Dog, Domestic")
SpeciesActiveModel("Turkey")
SpeciesActiveModel("Bear")
SpeciesActiveModel("Wolf")
SpeciesActiveModel("Opossum")
SpeciesActiveModel("Grouse")
SpeciesActiveModel("Porcupine")
SpeciesActiveModel("Fisher")
SpeciesActiveModel("Bobcat")
SpeciesActiveModel("Fox, Red")
SpeciesActiveModel("Mink")
SpeciesActiveModel("Otter")
SpeciesActiveModel("Cat, Domestic")
SpeciesActiveModel("Skunk, Striped")
SpeciesActiveModel("Woodchuck")
SpeciesActiveModel("Crane, Sandhill")
SpeciesActiveModel("Pheasant")
SpeciesActiveModel("Weasel")
SpeciesActiveModel("Badger")
SpeciesActiveModel("Elk")
SpeciesActiveModel("Beaver")
SpeciesActiveModel("Fox, Gray")
SpeciesActiveModel("Snowshoe Hare")
SpeciesActiveModel("Other Small Mammal")
SpeciesActiveModel("Muskrat")
SpeciesActiveModel("Other Domestic")
SpeciesActiveModel("Cougar")
SpeciesActiveModel("Moose")
SpeciesActiveModel("Marten")
SpeciesActiveModel("Crane, Whooping")
SpeciesActiveModel("Reptiles and Amphibians")
SpeciesActiveModel("Skunk, Spotted")
SpeciesActiveModel("Jackrabbit")

# How nocturnal is each species? -----------------------------------------------
# what % of detections occur at night?

sp <- unique(detections$SPECIES) # a list of species
Diel <- data_frame(species=character(), # an empty dataframe
                   nocturnality=double())
for(i in 1:length(sp)) {
  df <- subset(detections, detections$SPECIES==sp[i])
  x <- data_frame(species = sp[i])
  x$nocturnality <- nrow(subset(df, df$dielCategory=="NIGHT"))/nrow(df) # nocturnal detections/ total detections
  Diel <- rbind(Diel, x)
}

 # Is diel activity correlated with night light? --------------------------------

DielPercentActive <- 
  detections %>% 
  group_by(CAMERA_ID_CORRECTED_FOR_DUPLICATES, round(WI_Anthropocene_Mu_et_al_mNTL_Max, digits = 0), dielCategory, SPECIES) %>% # groups sites with similar NTL values together
  dplyr::summarise(count = n())  %>% # count the number of detections in each group
  pivot_wider(names_from = dielCategory, values_from = count)
DielPercentActive[is.na(DielPercentActive)] <- 0 # NAs should be zeros = sites with no detections
DielPercentActive$percentNocturnal <- DielPercentActive$NIGHT/(DielPercentActive$DAY + DielPercentActive$NIGHT)
names(DielPercentActive)[names(DielPercentActive) == 'round(WI_Anthropocene_Mu_et_al_mNTL_Max, digits = 0)'] <- 'NTL' #make that name prettier


NTL <- data_frame(species=character(),
                  NTL_rho=double(),
                  NTL_p=double())
for(i in 1:length(sp)) {
  df <- subset(DielPercentActive, DielPercentActive$SPECIES==sp[i])
  x <- data_frame(species = sp[i])
  x$NTL_rho <- cor.test(df$NTL, df$percentNocturnal, method = "spearman")$estimate # the correlation between NTL and nocturnality
  x$NTL_p <- cor.test(df$NTL, df$percentNocturnal, method = "spearman")$p.value # the significance of that correlation
  NTL <- rbind(NTL, x)
}
Diel <- merge(Diel, NTL, all=TRUE)

#plotting by species

NTLNocModel <- function(species) {
  species_diel_data = subset(DielPercentActive, SPECIES == species)
  ggplot(species_diel_data, aes(x = NTL, y = percentNocturnal)) +
    geom_point() +
    geom_smooth() +
    ggtitle(paste(species))
  ggsave(paste0(species, "_NTLNocModel.png"), path = 'C:/Users/amajerus/Desktop/Annabelle/Figures/NTLNocModel/')
  
}

NTLNocModel("Raccoon")
NTLNocModel("Coyote")
NTLNocModel("Deer")
NTLNocModel("Cottontail")
NTLNocModel("Other Bird")
NTLNocModel("Squirrels and Chipmunks")
NTLNocModel("Dog, Domestic")
NTLNocModel("Raccoon")
NTLNocModel("Turkey")
NTLNocModel("Bear")
NTLNocModel("Wolf")
NTLNocModel("Opossum")
NTLNocModel("Grouse")
NTLNocModel("Porcupine")
NTLNocModel("Fisher")
NTLNocModel("Bobcat")
NTLNocModel("Fox, Red")
NTLNocModel("Mink")
NTLNocModel("Otter")
NTLNocModel("Cat, Domestic")
NTLNocModel("Skunk, Striped")
NTLNocModel("Woodchuck")
NTLNocModel("Crane, Sandhill")
NTLNocModel("Pheasant")
NTLNocModel("Weasel")
NTLNocModel("Badger")
NTLNocModel("Elk")
NTLNocModel("Beaver")
NTLNocModel("Fox, Gray")
NTLNocModel("Snowshoe Hare")
NTLNocModel("Other Small Mammal")
NTLNocModel("Muskrat")
NTLNocModel("Other Domestic")
NTLNocModel("Cougar")
NTLNocModel("Moose")
NTLNocModel("Marten")
NTLNocModel("Crane, Whooping")
NTLNocModel("Reptiles and Amphibians")
NTLNocModel("Skunk, Spotted")
NTLNocModel("Jackrabbit")

# Is diel activity correlated with GHM? ----------------------------------------

DielGHMActive <- 
  detections %>% 
  group_by(round(WI_Anthropocene_NASA_SEDAC_nsGHM, 2), dielCategory, SPECIES) %>% 
  dplyr::summarise(count = n())  %>%
  pivot_wider(names_from = dielCategory, values_from = count)
DielGHMActive[is.na(DielGHMActive)] <- 0
DielGHMActive$percentNocturnal <- DielGHMActive$NIGHT/(DielGHMActive$DAY + DielGHMActive$NIGHT)
names(DielGHMActive)[names(DielGHMActive) == 'round(WI_Anthropocene_NASA_SEDAC_nsGHM, 2)'] <- 'GHM'

GHM <- data_frame(species=character(),
                  GHM_rho=double(),
                  GHM_p=double())
for(i in 1:length(sp)) {
  df <- subset(DielGHMActive, DielGHMActive$SPECIES==sp[i])
  x <- data_frame(species = sp[i])
  x$GHM_rho <- cor.test(df$GHM, df$percentNocturnal, method = "spearman")$estimate
  x$GHM_p <- cor.test(df$GHM, df$percentNocturnal, method = "spearman")$p.value
  GHM <- rbind(GHM, x)
}
Diel <- merge(Diel, GHM, all=TRUE)

GHMNocModel = function(species) {
  species_diel_data = subset(DielGHMActive, SPECIES == species)
  ggplot(species_diel_data, aes(x=GHM, y=percentNocturnal)) +
    geom_point() +
    geom_smooth() +
    ggtitle(print(species))
  ggsave(paste0(species, "_GHMNocModel.png"), path = 'C:/Users/amajerus/Desktop/Annabelle/Figures/GHMNocModel/')
  
}

GHMNocModel("Raccoon")
GHMNocModel("Coyote")
GHMNocModel("Deer")
GHMNocModel("Cottontail")
GHMNocModel("Other Bird")
GHMNocModel("Squirrels and Chipmunks")
GHMNocModel("Dog, Domestic")
GHMNocModel("Raccoon")
GHMNocModel("Turkey")
GHMNocModel("Bear")
GHMNocModel("Wolf")
GHMNocModel("Opossum")
GHMNocModel("Grouse")
GHMNocModel("Porcupine")
GHMNocModel("Fisher")
GHMNocModel("Bobcat")
GHMNocModel("Fox, Red")
GHMNocModel("Mink")
GHMNocModel("Otter")
GHMNocModel("Cat, Domestic")
GHMNocModel("Skunk, Striped")
GHMNocModel("Woodchuck")
GHMNocModel("Crane, Sandhill")
GHMNocModel("Pheasant")
GHMNocModel("Weasel")
GHMNocModel("Badger")
GHMNocModel("Elk")
GHMNocModel("Beaver")
GHMNocModel("Fox, Gray")
GHMNocModel("Snowshoe Hare")
GHMNocModel("Other Small Mammal")
GHMNocModel("Muskrat")
GHMNocModel("Other Domestic")
GHMNocModel("Cougar")
GHMNocModel("Moose")
GHMNocModel("Marten")
GHMNocModel("Crane, Whooping")
GHMNocModel("Reptiles and Amphibians")
GHMNocModel("Skunk, Spotted")
GHMNocModel("Jackrabbit")

#Is diel activity correlated with population density?

DielPopDActive <- 
  detections %>% 
  group_by(round(WI_Anthropocene_GPW_populationDensity, 2), dielCategory, SPECIES) %>% 
  dplyr::summarise(count = n())  %>%
  pivot_wider(names_from = dielCategory, values_from = count)
DielPopDActive[is.na(DielPopDActive)] <- 0
DielPopDActive$percentNocturnal <- DielPopDActive$NIGHT/(DielPopDActive$DAY + DielPopDActive$NIGHT)
names(DielPopDActive)[names(DielPopDActive) == 'round(WI_Anthropocene_GPW_populationDensity, 2)'] <- 'PopDActive'

PopDActive <- data_frame(species=character(),
                  PopDActive_rho=double(),
                  PopDActive_p=double())
for(i in 1:length(sp)) {
  df <- subset(DielPopDActive, DielPopDActive$SPECIES==sp[i])
  x <- data_frame(species = sp[i])
  x$PopDActive_rho <- cor.test(df$PopDActive, df$percentNocturnal, method = "spearman")$estimate
  x$PopDActive_p <- cor.test(df$PopDActive, df$percentNocturnal, method = "spearman")$p.value
  PopDActive <- rbind(PopDActive, x)
}
Diel <- merge(Diel, PopDActive, all=TRUE)

PopDActiveNocModel = function(species) {
  species_diel_data = subset(DielPopDActive, SPECIES == species)
  ggplot(species_diel_data, aes(x=PopDActive, y=percentNocturnal)) +
    geom_point() +
    geom_smooth() +
    ggtitle(print(species))
  ggsave(paste0(species, "_PopDactiveNocModel.png"), path = 'C:/Users/amajerus/Desktop/Annabelle/Figures/PopDActiveNocModel/')
  
}

PopDActiveNocModel("Raccoon")
PopDActiveNocModel("Coyote")
PopDActiveNocModel("Deer")
PopDActiveNocModel("Cottontail")
PopDActiveNocModel("Other Bird")
PopDActiveNocModel("Squirrels and Chipmunks")
PopDActiveNocModel("Dog, Domestic")
PopDActiveNocModel("Raccoon")
PopDActiveNocModel("Turkey")
PopDActiveNocModel("Bear")
PopDActiveNocModel("Wolf")
PopDActiveNocModel("Opossum")
PopDActiveNocModel("Grouse")
PopDActiveNocModel("Porcupine")
PopDActiveNocModel("Fisher")
PopDActiveNocModel("Bobcat")
PopDActiveNocModel("Fox, Red")
PopDActiveNocModel("Mink")
PopDActiveNocModel("Otter")
PopDActiveNocModel("Cat, Domestic")
PopDActiveNocModel("Skunk, Striped")
PopDActiveNocModel("Woodchuck")
PopDActiveNocModel("Crane, Sandhill")
PopDActiveNocModel("Pheasant")
PopDActiveNocModel("Weasel")
PopDActiveNocModel("Badger")
PopDActiveNocModel("Elk")
PopDActiveNocModel("Beaver")
PopDActiveNocModel("Fox, Gray")
PopDActiveNocModel("Snowshoe Hare")
PopDActiveNocModel("Other Small Mammal")
PopDActiveNocModel("Muskrat")
PopDActiveNocModel("Other Domestic")
PopDActiveNocModel("Cougar")
PopDActiveNocModel("Moose")
PopDActiveNocModel("Marten")
PopDActiveNocModel("Crane, Whooping")
PopDActiveNocModel("Reptiles and Amphibians")
PopDActiveNocModel("Skunk, Spotted")
PopDActiveNocModel("Jackrabbit")

#Building a model for access to cities

DielCityActive <- 
  detections %>% 
  group_by(round(WI_Anthropocene_Weiss_et_al_A2C, 2), dielCategory, SPECIES) %>% 
  dplyr::summarise(count = n())  %>%
  pivot_wider(names_from = dielCategory, values_from = count)
DielCityActive[is.na(DielCityActive)] <- 0
DielCityActive$percentNocturnal <- DielCityActive$NIGHT/(DielCityActive$DAY + DielCityActive$NIGHT)
names(DielCityActive)[names(DielCityActive) == 'round(WI_Anthropocene_Weiss_et_al_A2C, 2)'] <- 'City'

City <- data_frame(species=character(),
                  City_rho=double(),
                  City_p=double())
for(i in 1:length(sp)) {
  df <- subset(DielCityActive, DielCityActive$SPECIES==sp[i])
  x <- data_frame(species = sp[i])
  x$City_rho <- cor.test(df$City, df$percentNocturnal, method = "spearman")$estimate
  x$City_p <- cor.test(df$City, df$percentNocturnal, method = "spearman")$p.value
  City <- rbind(City, x)
}
Diel <- merge(Diel, City, all=TRUE)

CityNocModel = function(species) {
  species_diel_data = subset(DielCityActive, SPECIES == species)
  ggplot(species_diel_data, aes(x=City, y=percentNocturnal)) +
    geom_point() +
    geom_smooth() +
    ggtitle(print(species))
  ggsave(paste0(species, "_CityNocModel.png"), path = 'C:/Users/amajerus/Desktop/Annabelle/Figures/CityNocModel/')
  
}

CityNocModel("Raccoon")
CityNocModel("Coyote")
CityNocModel("Deer")
CityNocModel("Cottontail")
CityNocModel("Other Bird")
CityNocModel("Squirrels and Chipmunks")
CityNocModel("Dog, Domestic")
CityNocModel("Raccoon")
CityNocModel("Turkey")
CityNocModel("Bear")
CityNocModel("Wolf")
CityNocModel("Opossum")
CityNocModel("Grouse")
CityNocModel("Porcupine")
CityNocModel("Fisher")
CityNocModel("Bobcat")
CityNocModel("Fox, Red")
CityNocModel("Mink")
CityNocModel("Otter")
CityNocModel("Cat, Domestic")
CityNocModel("Skunk, Striped")
CityNocModel("Woodchuck")
CityNocModel("Crane, Sandhill")
CityNocModel("Pheasant")
CityNocModel("Weasel")
CityNocModel("Badger")
CityNocModel("Elk")
CityNocModel("Beaver")
CityNocModel("Fox, Gray")
CityNocModel("Snowshoe Hare")
CityNocModel("Other Small Mammal")
CityNocModel("Muskrat")
CityNocModel("Other Domestic")
CityNocModel("Cougar")
CityNocModel("Moose")
CityNocModel("Marten")
CityNocModel("Crane, Whooping")
CityNocModel("Reptiles and Amphibians")
CityNocModel("Skunk, Spotted")
CityNocModel("Jackrabbit")

#Building model for impacting sound

DielSoundActive <- 
  detections %>% 
  group_by(round(WI_Anthropocene_NPS_impSound, 2), dielCategory, SPECIES) %>% 
  dplyr::summarise(count = n())  %>%
  pivot_wider(names_from = dielCategory, values_from = count)
DielSoundActive[is.na(DielSoundActive)] <- 0
DielSoundActive$percentNocturnal <- DielSoundActive$NIGHT/(DielSoundActive$DAY + DielSoundActive$NIGHT)
names(DielSoundActive)[names(DielSoundActive) == 'round(WI_Anthropocene_NPS_impSound, 2)'] <- 'Sound'

Sound <- data_frame(species=character(),
                   Sound_rho=double(),
                   Sound_p=double())
for(i in 1:length(sp)) {
  df <- subset(DielSoundActive, DielSoundActive$SPECIES==sp[i])
  x <- data_frame(species = sp[i])
  x$Sound_rho <- cor.test(df$Sound, df$percentNocturnal, method = "spearman")$estimate
  x$Sound_p <- cor.test(df$Sound, df$percentNocturnal, method = "spearman")$p.value
  Sound <- rbind(Sound, x)
}
Diel <- merge(Diel, Sound, all=TRUE)

SoundNocModel = function(species) {
  species_diel_data = subset(DielSoundActive, SPECIES == species)
  ggplot(species_diel_data, aes(x=Sound, y=percentNocturnal)) +
    geom_point() +
    geom_smooth() +
    ggtitle(print(species))
  ggsave(paste0(species, "_SoundNocModel.png"), path = 'C:/Users/amajerus/Desktop/Annabelle/Figures/SoundNocModel/')
  
}

SoundNocModel("Raccoon")
SoundNocModel("Coyote")
SoundNocModel("Deer")
SoundNocModel("Cottontail")
SoundNocModel("Other Bird")
SoundNocModel("Squirrels and Chipmunks")
SoundNocModel("Dog, Domestic")
SoundNocModel("Raccoon")
SoundNocModel("Turkey")
SoundNocModel("Bear")
SoundNocModel("Wolf")
SoundNocModel("Opossum")
SoundNocModel("Grouse")
SoundNocModel("Porcupine")
SoundNocModel("Fisher")
SoundNocModel("Bobcat")
SoundNocModel("Fox, Red")
SoundNocModel("Mink")
SoundNocModel("Otter")
SoundNocModel("Cat, Domestic")
SoundNocModel("Skunk, Striped")
SoundNocModel("Woodchuck")
SoundNocModel("Crane, Sandhill")
SoundNocModel("Pheasant")
SoundNocModel("Weasel")
SoundNocModel("Badger")
SoundNocModel("Elk")
SoundNocModel("Beaver")
SoundNocModel("Fox, Gray")
SoundNocModel("Snowshoe Hare")
SoundNocModel("Other Small Mammal")
SoundNocModel("Muskrat")
SoundNocModel("Other Domestic")
SoundNocModel("Cougar")
SoundNocModel("Moose")
SoundNocModel("Marten")
SoundNocModel("Crane, Whooping")
SoundNocModel("Reptiles and Amphibians")
SoundNocModel("Skunk, Spotted")
SoundNocModel("Jackrabbit")

#saving the 'Diel' Data as a csv
write.csv(Diel, "C:/Users/amajerus/Desktop/Annabelle/CSVs/Diel.csv", row.names = FALSE)

#Putting appropriate data into data frame for colinearity test

dielColinear = select(detections, 
                      #SPECIES, dielCategory, 
                      WI_Anthropocene_GPW_populationDensity, 
                      WI_Anthropocene_Li_et_al_NTL_2013,
                      WI_Anthropocene_NASA_SEDAC_nsGHM,
                      WI_Anthropocene_NPS_impSound,
                      WI_Anthropocene_Weiss_et_al_A2C) %>% 
  rename('Population Density' = WI_Anthropocene_GPW_populationDensity,
         'NTL'= WI_Anthropocene_Li_et_al_NTL_2013,
         'GHM' = WI_Anthropocene_NASA_SEDAC_nsGHM,
         'Impacting Sound' = WI_Anthropocene_NPS_impSound,
         'Access to Cities' = WI_Anthropocene_Weiss_et_al_A2C)


#inverting access to cities data
dielColinear$"Access to Cities" <- 1 / dielColinear$"Access to Cities"

library(corrplot)
correlation = cor(dielColinear, method = "spearman")

png(filename = "C:/Users/amajerus/Desktop/Annabelle/Figures/Figure_3.png", 
    width = 8.91, height = 5.74, units = "in", res = 100 )
corrplot(correlation,method = "color",
         cex.main = 0.8,
         mar=c(0,0,2,0),
         type = "upper", 
         tl.cex = 0.8,
         tl.col = "black",
  tl.srt = 30,
  addCoef.col = "black")
dev.off()



#making a facet plot
GHM_species = c("Wolf", "Raccoon", "Fisher", "Coyote", "Cottontail", "Squirrels and Chipmunks", "Porcupine", "Otter")
##getting the spearman values for anthropogenic factors of significant species
Diel_sig = Diel %>% select(-contains("_p")) %>% select(-contains("lunillum")) %>% select(-contains("UrbanActive")) %>% 
  rename("Impacting Sound" = "Sound_rho", "Population Density" = "PopDActive_rho", 
           "NTL" = "NTL_rho","GHM" = "GHM_rho","Access to Cities" = "City_rho") %>% 
  filter(species %in% GHM_species) %>% 
  pivot_longer(cols = c("GHM", "Population Density", "NTL", "GHM", "Access to Cities", "Impacting Sound"),
             names_to = "variable",
             values_to = "correlation")
#inverting access to cities
Diel_sig$correlation[Diel_sig$variable == "Access to Cities"] <- 
  -Diel_sig$correlation[Diel_sig$variable == "Access to Cities"]

# Define significant variables for each species
significant_vars <- list(
  "Squirrels and Chipmunks" = c("GHM", "Impacting Sound"),
  "Porcupine" = c("GHM"),
  "Otter" = c("GHM", "Access to Cities")
)

#Filtering based on significance
Diel_sig <- Diel_sig %>%
  rowwise() %>%
  filter(
    (species %in% names(significant_vars) && variable %in% significant_vars[[species]]) ||
      !(species %in% names(significant_vars))
  ) %>%
  ungroup()

#Setting Species as a factor to group together diet types
Diel_sig$species <- factor(Diel_sig$species, levels = c("Porcupine", 
                                                        "Squirrels and Chipmunks",
                                                        "Cottontail",
                                                        "Otter",
                                                        "Fisher",
                                                        "Raccoon",
                                                        "Coyote",
                                                        "Wolf"))
# Drop unused levels from the 'variable' factor
Diel_sig$variable <- factor(Diel_sig$variable, levels = unique(Diel_sig$variable))

#plotting
Spearman_sig = 
  ggplot(Diel_sig, aes(x=correlation, y=variable, fill = species)) +
  geom_col() +
  facet_grid(rows = vars(species)) +
  theme_minimal() +
  facet_grid(rows = vars(species), scales = "free_y", space = "free") +
  scale_fill_manual(values = c("Wolf" = "#DA6C6C",
                               "Raccoon" = "#7B4019",
                               "Fisher" = "#FF7D29",
                               "Coyote" = "#FFBF78",
                               "Cottontail" = "#a6d854",
                               "Squirrels and Chipmunks" = "#06923E",
                               "Porcupine" = "#5D8736",
                               "Otter" = "#EA2F14")) +       # pick your favorite colors https://colorhunt.co/
  geom_vline(xintercept=0) +                                        # add a line at 0
  xlim(-1, 1) +                                                     # set axis limits
  ylab("Geospatial Variable") +                                  # custom axes labels
  xlab("Correlation Coefficient")
Spearman_sig
ggsave("Spearman_sig.png", path = "C:\\Users\\amajerus\\Desktop\\Annabelle\\Figures\\Spearmans")

#plotting ghm vs percent nocturnal for significant species on the same plot

# Filter the dataset for the significant species
GHM_sig <- subset(DielGHMActive, SPECIES %in% GHM_species)
GHM_sig$SPECIES <- factor(GHM_sig$SPECIES, levels = c("Porcupine", 
                                                        "Squirrels and Chipmunks",
                                                        "Cottontail",
                                                        "Otter",
                                                        "Fisher",
                                                        "Raccoon",
                                                        "Coyote",
                                                        "Wolf"))
ghm_color = c("Wolf" = "#DA6C6C",
    "Raccoon" = "#7B4019",
    "Fisher" = "#FF7D29",
    "Coyote" = "#FFBF78",
    "Cottontail" = "#a6d854",
    "Squirrels and Chipmunks" = "#06923E",
    "Porcupine" = "#5D8736",
    "Otter" = "#EA2F14")


# Create the scatter plot with trend lines
GHMNocModel_ScatterPlot = ggplot(GHM_sig, aes(x = GHM, y = percentNocturnal, color = SPECIES)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", se = FALSE) +  # You can change method to "loess" if preferred
  labs(title = "GHM vs Percent Nocturnal by Species",
       x = "GHM",
       y = "Proportion Nocturnal",
       color = "Species") +
  scale_color_manual(values = ghm_color) +
  theme_minimal() 
GHMNocModel_ScatterPlot
ggsave("GHMNocModel_ScatterPlot.png", path = "C:/Users/amajerus/Desktop/Annabelle/Figures/GHMNocModel/")

#Creating a map of camera density for snapshot wisconsin project
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)

wi_counties = st_read("C:\\Users\\amajerus\\Desktop\\Annabelle\\CSVs\\County_Boundaries_24K.geojson")
cameras = read.csv("C:\\Users\\amajerus\\Desktop\\Annabelle\\CSVs\\camera_ID_corrected_for_duplicates.csv")

camera_density = ggplot(cameras) +
  geom_point(data = cameras, aes(x = LON, y = LAT), color = "#E83F25", size = 1.25, alpha = 0.4) +
  geom_sf(data = wi_counties, fill = "transparent", color = "black", linewidth = 0.6) + 
  theme_light() + 
  ylab("Latitude") +
  xlab("Longitude")
camera_density
ggsave("Camera_density.png", path = "C:/Users/amajerus/Desktop/Annabelle/Figures/")

#total observations csv

observation_df = detections %>% group_by(SPECIES) %>% summarise(number_obs = n())
write.csv(observation_df,"C:\\Users\\amajerus\\Desktop\\Annabelle\\CSVs\\observations.csv", row.names = FALSE )

#plotting biocube variables on map
biocube_variables = detections[,c("COUNTY_NAME.y", "WI_Anthropocene_GPW_populationDensity", 
                                  "WI_Anthropocene_Li_et_al_NTL_2013", "WI_Anthropocene_NASA_SEDAC_nsGHM",
                                  "WI_Anthropocene_NPS_impSound",
                                  "WI_Anthropocene_Weiss_et_al_A2C", "LON", "LAT")] %>% rename("COUNTY_NAME" = "COUNTY_NAME.y") %>% 
  distinct(COUNTY_NAME, .keep_all = TRUE) 

biocube_variables$COUNTY_NAME = str_replace(biocube_variables$COUNTY_NAME, "Fond Du Lac", "Fond du Lac")

biowi_map_data = wi_counties %>% left_join(biocube_variables, by = "COUNTY_NAME")

#ghm
ghm_map = ggplot(biowi_map_data) +
  geom_sf(aes(fill = WI_Anthropocene_NASA_SEDAC_nsGHM)) +
  geom_sf_text(aes(label = round(WI_Anthropocene_NASA_SEDAC_nsGHM, 2)), size = 2, color = "black", fontface = "bold") +
  scale_fill_gradient(low = "#FFAF45", high = "#211C84", na.value = "grey90", 
                      limits = c(0, 1),
                      breaks = c(0, 0.25, 0.5, 0.75, 1),
                      labels = c("0", "0.25", "0.5", "0.75", "1")
  ) +
  labs(
       fill = "Cumulative GHM") +
  theme_minimal() + xlab("Longitude") +ylab("Latitude")
ghm_map
ggsave("ghm_map.png", path = "C:\\Users\\amajerus\\Desktop\\Annabelle\\Figures\\Maps")

#popD
# Calculate centroids for placing bubbles and labels
biowi_centroids <- biowi_map_data %>%
  st_centroid() %>%
  st_coordinates() %>%
  as.data.frame() %>%
  bind_cols(biowi_map_data)

# Plot with bubbles and labels
popD_map <- ggplot(biowi_map_data) +
  geom_sf(fill = "white", color = "black") +  # base map
  geom_point(data = biowi_centroids,
             aes(x = X, y = Y, size = WI_Anthropocene_GPW_populationDensity),
             color = "#FF7601", alpha = 0.6) +
  geom_text(data = biowi_centroids,
            aes(x = X, y = Y, label = round(WI_Anthropocene_GPW_populationDensity, 1)),
            size = 2.5, color = "black", vjust = 0) +  # adjust vjust to position labels
  scale_size_continuous(range = c(1, 12), limits = c(0, 1500),
                        breaks = c(0, 500, 1000, 1500),
                        labels = c("0", "500", "1000", "1500")) +
  labs(
       size = bquote("Number of People per Mile"^2)) +
  theme_minimal()+ xlab("Longitude") +ylab("Latitude")
popD_map
ggsave("popD_map.png", path = "C:\\Users\\amajerus\\Desktop\\Annabelle\\Figures\\Maps")


#NTL
NTL_map = ggplot(biowi_map_data) +
  geom_sf(aes(fill = WI_Anthropocene_Li_et_al_NTL_2013)) +
  geom_sf_text(aes(label = round(WI_Anthropocene_Li_et_al_NTL_2013, 2)), size = 2, color = "black", 
               fontface = "bold"
  ) +
  scale_fill_gradient(low = "#9FC87E", high = "#8A0000", na.value = "grey80", 
                      limits = c(0, 63),
                      breaks = c(0, 31.5, 63),
                      labels = c("0", "31.5","63")
  ) +
  labs(,
       fill = "NTL (DN Value)") + xlab("Longitude") +ylab("Latitude") +
  theme_minimal()
NTL_map
ggsave("NTL_map.png", path = "C:\\Users\\amajerus\\Desktop\\Annabelle\\Figures\\Maps")

#impacting sound
NTL_map = ggplot(biowi_map_data) +
  geom_sf(aes(fill = WI_Anthropocene_NPS_impSound)) +
  geom_sf_text(aes(label = round(WI_Anthropocene_NPS_impSound, 2)), size = 2, color = "black", 
               fontface = "bold"
  ) +
  scale_fill_gradient(low = "#77BEF0", high = "#FE5D26", na.value = "grey80", 
                      limits = c(0, 20),
                      breaks = c(0, 5, 10, 15, 20),
                      labels = c("0", "5","10", "15", "20")
  ) +
  labs(
       fill = "Impacting Sound (L50, dBA)") + xlab("Longitude") +ylab("Latitude") +
  theme_minimal()
NTL_map
ggsave("Sound_map.png", path = "C:\\Users\\amajerus\\Desktop\\Annabelle\\Figures\\Maps")

#Access to cities
city_map = ggplot(biowi_map_data) +
  geom_sf(aes(fill = WI_Anthropocene_Weiss_et_al_A2C)) +
  geom_sf_text(aes(label = round(WI_Anthropocene_Weiss_et_al_A2C, 2)), size = 2, color = "black", 
               fontface = "bold"
  ) +
  scale_fill_gradient(low = "yellow", high = "#F72C5B", na.value = "grey80", 
                      limits = c(0, 200),
                      breaks = c(0, 50, 100, 150, 200),
                      labels = c("0", "50","100", "150", "200")
  ) +
  labs(
       fill = "Travel Time to City (min)") + xlab("Longitude") +ylab("Latitude") +
  theme_minimal()
city_map
ggsave("city_map.png", path = "C:\\Users\\amajerus\\Desktop\\Annabelle\\Figures\\Maps")

