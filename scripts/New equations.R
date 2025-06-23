### Define new equations #####

library(isocat)
library(terra)
library(tidyverse)
library(mapview)
library(readr)
library(rnaturalearth)
library(RColorBrewer)
library(rnaturalearth)

## prepare environment
rm(list=ls()) # clear the global environment to save on space
# rm(object_name) # to remove specific objects
set.seed(725)
select <- dplyr::select
custom.palette <- colorRampPalette(brewer.pal(9, "YlGnBu")) # Color palette
gc()

# Need shapefiles for WI and MN to clip out the samples included in the PLOS dataset
US.states <- ne_states(country = "United States of America", returnclass = "sv") 
WIMN <- US.states[US.states$name %in% c("Wisconsin","Minnesota"),]

## Load in data needed to run the spatial assignment
# dataframe with sample d2h values and metadata
new.knownorig <- read_rds("data/processed/cleaned_feather_results/reference_data_df.rds") %>% 
  filter(age == "local") %>% # for now, not using AHY or HY birds 
  filter(!(state == "WI" & species == "RNDU")) %>% # the two rndu samples from WI fall outside eBird breeding range crop
  mutate(species = tolower(species)) %>%
  rename(dataset = "state",
         d2H = "d2h_vs_vsmow_original") %>%
  mutate(species = toupper(species)) %>%
  select(dataset, species, longitude, latitude, d2H) %>%
  mutate(ifelse(species == "RNDU", 'Divers', 'Dabblers')) %>%
  rename(guild = 'ifelse(species == "RNDU", "Divers", "Dabblers")') %>%
  vect(geom = c("longitude","latitude"), crs = "EPSG:4326")
  
knownorig <- read.csv("data/raw/kusack_plos_2023/kusack_2023_known_origin_datasets_d2H.csv") %>%
  vect(geom = c("long","lat"), crs = "EPSG:4326") 
knownorig <- knownorig[!(terra::is.related(knownorig, buffer(WIMN, 100), relation = "coveredby"))] # need the buffer to catch border samples
  
knownorig <- rbind(knownorig, new.knownorig)

mapview(knownorig)

# growing-season isoscape
isoscape_gs <- rast("data/raw/precipitation_isoscapes/gs_isoscape/d2h_GS.tif")
crs(isoscape_gs) <- "EPSG:4326"

# mean annual isoscape
isoscape_ma <- rast("data/raw/precipitation_isoscapes/ma_isoscape/d2h_MA.tif")
crs(isoscape_ma) <- "EPSG:4326"




knownorig$gsd <- terra::extract(isoscape_gs, knownorig)[,2]
knownorig$mad <- terra::extract(isoscape_ma, knownorig)[,2]

## GSD

eq.df <- data.frame(guild = c("Dabblers","Divers"),  equation = c(NA,NA))
  
eq.df$equation[eq.df$guild == "Dabblers"] <- paste0("delta^2*H[f]~`=`~", 
                                                   round(coef(lm(d2H ~ gsd, data = knownorig[knownorig$guild == "Dabblers",]))[1],2), "~+~", 
                                                   round(coef(lm(d2H ~ gsd, data = knownorig[knownorig$guild == "Dabblers",]))[2],2), "~x~delta^2*H[p]")

eq.df$equation[eq.df$guild == "Divers"] <- paste0("delta^2*H[f]~`=`~", 
                                                   round(coef(lm(d2H ~ gsd, data = knownorig[knownorig$guild == "Divers",]))[1],2), "~+~", 
                                                   round(coef(lm(d2H ~ gsd, data = knownorig[knownorig$guild == "Divers",]))[2],2), "~x~delta^2*H[p]")

(p.gsd <- ggplot(data.frame(knownorig), aes(x = gsd, y = d2H)) + 
  geom_point(aes(col = dataset)) + 
  facet_wrap(~guild) + 
  stat_smooth(method = "lm") +
  geom_text(data = eq.df, x = -80, y = -250, aes(label = equation), check_overlap = TRUE,
            fontface = "italic", parse = T, size = 2) + 
  labs(col = "Dataset", y = "Feather d2H", x = "Precip. d2H (gsd)") + 
  theme_classic())

ggsave(filename = "figures/gsd.equations.png", p.gsd, dpi = 300, width = 6, height = 4.5, units = "in")

## MAD

eq.df <- data.frame(guild = c("Dabblers","Divers"),  equation = c(NA,NA))

eq.df$equation[eq.df$guild == "Dabblers"] <- paste0("delta^2*H[f]~`=`~", 
                                                    round(coef(lm(d2H ~ mad, data = knownorig[knownorig$guild == "Dabblers",]))[1],2), "~+~", 
                                                    round(coef(lm(d2H ~ mad, data = knownorig[knownorig$guild == "Dabblers",]))[2],2), "~x~delta^2*H[p]")

eq.df$equation[eq.df$guild == "Divers"] <- paste0("delta^2*H[f]~`=`~", 
                                                  round(coef(lm(d2H ~ mad, data = knownorig[knownorig$guild == "Divers",]))[1],2), "~+~", 
                                                  round(coef(lm(d2H ~ mad, data = knownorig[knownorig$guild == "Divers",]))[2],2), "~x~delta^2*H[p]")

(p.mad <- ggplot(data.frame(knownorig), aes(x = mad, y = d2H)) + 
  geom_point(aes(col = dataset)) + 
  facet_wrap(~guild) + 
  stat_smooth(method = "lm") +
  geom_text(data = eq.df, x = -80, y = -250, aes(label = equation), check_overlap = TRUE,
            fontface = "italic", parse = T, size = 2) + 
  labs(col = "Dataset", y = "Feather d2H", x = "Precip. d2H (mad)") + 
  theme_classic())

ggsave(filename = "figures/mad.equations.png", p.mad, dpi = 300, width = 6, height = 4.5, units = "in")




