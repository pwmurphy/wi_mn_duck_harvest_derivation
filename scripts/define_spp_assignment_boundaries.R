######################## defining origin assignment boundaries for MALL, RNDU, WODU ################################
######### using spp range maps from eBird and BirdLife International clipped to the Miss and Central flyways #######

# libraries
library(sf)
library(rmapshaper)
library(ebirdst)
library(rnaturalearth)
library(tidyverse)
library(janitor)
library(mapview)
library(fs)

select <- dplyr::select

### isolate the mississippi and central flyways to clip range boundaries to ###

# load north america
northamerica <- ne_countries(continent = "North America", scale = 50) %>% 
  filter(!admin == "Greenland") # remove Greenland
mapview(northamerica)

# Canada
# partition Canadian provinces to enable manually formation of respective flyways
canada <- ne_states(country= "canada")
manitoba <- canada %>% filter(name_en == "Manitoba")
ontario <- canada %>% filter(name_en == "Ontario")
nunavut <- canada %>% filter(name_en == "Nunavut")
saskatchewan <- canada %>% filter(name_en == "Saskatchewan")
alberta <- canada %>% filter(name_en == "Alberta")
north_terr <- canada %>% filter(name_en == "Northwest Territories")

# bind individual provinces into respective flyways
ms_flyway_canada <- rbind(manitoba,ontario,nunavut) %>% 
  select(name)
mapview(ms_flyway_canada)

cn_flyway_canada<-rbind(saskatchewan, alberta, north_terr) %>% 
  select(name)
mapview(cn_flyway_canada)

# US
# Read in all four waterfowl flyways 
us_flyways <- st_read("data/raw/flyways/us_flyways.shp") %>% 
  st_transform(4326) %>% 
  ms_simplify(keep=0.05) %>% 
  clean_names()

mapview(us_flyways)

# Subset all US flyways down to independent flyways
ms_us_flyway <- us_flyways %>% filter(name == "Mississippi Flyway")
cn_us_flyway <- us_flyways %>% filter(name == "Central Flyway")

# merge US flyway shapefiles with respective Canadian flyways
ms_contin_flyway <- ms_us_flyway %>% 
  st_transform(st_crs(ms_flyway_canada)) %>% 
  rbind(ms_flyway_canada) %>% 
  mutate(name = case_when(
    name == "Mississippi Flyway" ~ "United States",  # Change "Mississippi" to "United States"
    TRUE ~ name  # Keep all other values the same
  ))
mapview(ms_contin_flyway)

cn_contin_flyway <- cn_us_flyway %>% 
  st_transform(st_crs(cn_flyway_canada)) %>% 
  rbind(cn_flyway_canada) %>% 
  mutate(name = case_when(
    name == "Central Flyway" ~ "United States",  # Change "Mississippi" to "United States"
    TRUE ~ name  # Keep all other values the same
  ))

mapview(cn_contin_flyway)

# save flyways
write_rds(ms_contin_flyway, file = "data/processed/isolated_flyways/ms_contin_flyway.rds")
write_rds(cn_contin_flyway, file = "data/processed/isolated_flyways/cn_contin_flyway.rds")

# merge continental flyways together
ms_cn_flyway <- rbind(ms_contin_flyway, cn_contin_flyway)
mapview(ms_cn_flyway)
# save
write_rds(ms_cn_flyway, "data/processed/isolated_flyways/ms_cn_contin_flyway.rds")

### download and crop eBird ranges ###

# # set eBird key acquired through personal ebird account (only need to set onces)
# set_ebirdst_access_key("792o0s8afe6k")

# # initial data download (only do once, takes several minutes)

# # get species name from database
# mallard <- ebirdst_runs %>% filter(common_name == "Mallard")
# rndu <- ebirdst_runs %>% filter(common_name == "Ring-necked Duck")
# wodu <- ebirdst_runs %>% filter(common_name == "Wood Duck")
# 
# # download data
# ebirdst_download_status(species = "mallar3", download_ranges = TRUE)
# ebirdst_download_status(species = "rinduc", download_ranges = TRUE)
# ebirdst_download_status(species = "wooduc", download_ranges = TRUE)
# 
# # where ebird data is stored
# ebird <- ebirdst_data_dir()

# list all the files I have downloaded from ebird
dir_ls(ebirdst_data_dir(), recurse = TRUE) %>%
  as_tibble() %>%
  print(n = Inf)

## load in the mall, rndu, wodu breeding ranges 
mall_range_ebird <- load_ranges(species = "mallar3", resolution = "27km", smoothed = TRUE, path = ebirdst_data_dir()) %>% filter(season == 'breeding')
wodu_range_ebird <- load_ranges(species = "wooduc", resolution = "27km", smoothed = TRUE, path = ebirdst_data_dir()) %>% filter(season == 'breeding')
rndu_range_ebird <- load_ranges(species = "rinduc", resolution = "27km", smoothed = TRUE, path = ebirdst_data_dir()) %>% filter(season == 'breeding')

# save the full global ranges
write_rds(mall_range_ebird, "data/processed/spp_ranges/mall_ebird_range.rds")
write_rds(wodu_range_ebird, "data/processed/spp_ranges/wodu_ebird_range.rds")
write_rds(rndu_range_ebird, "data/processed/spp_ranges/rndu_ebird_range.rds")

## trim global ranges to Mississippi and Central Flyways

# Restrict breeding range to North America (minus Greenland)
mall_range_ebird_NA <- st_intersection(mall_range_ebird, northamerica) %>% summarise()
wodu_range_ebird_NA <- st_intersection(wodu_range_ebird, northamerica) %>% summarise()
rndu_range_ebird_NA <- st_intersection(rndu_range_ebird, northamerica) %>% summarise()

# Clip breeding ranges down to the Miss and Central Flyway
mall_ebird_assignment_range <- st_intersection(mall_range_ebird_NA, ms_cn_flyway) %>%
  summarise()
mapview(mall_ebird_assignment_range)

wodu_ebird_assignment_range <- st_intersection(wodu_range_ebird_NA, ms_cn_flyway) %>%
  summarise()
mapview(wodu_ebird_assignment_range)

rndu_ebird_assignment_range <- st_intersection(rndu_range_ebird_NA, ms_cn_flyway) %>%
  summarise()
mapview(rndu_ebird_assignment_range)

# save assignment ranges
write_rds(mall_ebird_assignment_range, "data/processed/spp_assignment_ranges/mall_ebird_assignment_range.rds")
write_rds(wodu_ebird_assignment_range, "data/processed/spp_assignment_ranges/wodu_ebird_assignment_range.rds")
write_rds(rndu_ebird_assignment_range, "data/processed/spp_assignment_ranges/rndu_ebird_assignment_range.rds")

###########
### download and crop BirdLife International ranges ###

## requested range data from the BLI website-- had to fill out short survey, data download is a zip file with range maps for all species (17k +) so it's a big file
# data came as a gpkg with a Word doc for metadata that gives fields and value explanations

# big file to load in at once-- instead, subset for species we want
gpkg_file <- "C:/Users/murphpwx/Documents/bli_ranges_all_spp/BOTW_2024_2.gpkg"
st_layers(gpkg_file) # look at the layers in it, only want the 'all_species' layer 

# # want to subset to our species, pull sample values to see how they are formatted 
# sci_name_preview <- st_read(gpkg_file, layer = "all_species", 
#                             query = "SELECT DISTINCT sci_name FROM all_species LIMIT 10")
# print(sci_name_preview)
# species names in sci_name variable and formatted like "Genus species"
# range type (breeding vs nonbreeding vs year round stored in seasonal variable, 2 = breeding range and resident = 1)

# subset the dataset as it's being read with a SQL query

mall_range_bli <- st_read(gpkg_file, 
                          query = "SELECT * FROM all_species 
                                       WHERE sci_name = 'Anas platyrhynchos' 
                                       AND seasonal IN (1,2)") %>% 
  st_make_valid() # mallard range has invalid geometry 
mapview(mall_range_bli)

wodu_range_bli <- st_read(gpkg_file, 
                      query = "SELECT * FROM all_species 
                                       WHERE sci_name = 'Aix sponsa' 
                                       AND seasonal IN (1,2)")
mapview(wodu_range_bli)

rndu_range_bli <- st_read(gpkg_file, 
                      query = "SELECT * FROM all_species 
                                       WHERE sci_name = 'Aythya collaris' 
                                       AND seasonal IN (1,2)")
mapview(rndu_range_bli)

# save global ranges 
write_rds(mall_range_bli, "data/processed/spp_ranges/mall_bli_range.rds")
write_rds(wodu_range_bli, "data/processed/spp_ranges/wodu_bli_range.rds")
write_rds(rndu_range_bli, "data/processed/spp_ranges/rndu_bli_range.rds")

## trim global ranges to Mississippi and Central Flyways

# Restrict breeding range to North America (minus Greenland)
mall_range_bli_NA <- st_intersection(mall_range_bli, northamerica) %>% summarise()
wodu_range_bli_NA <- st_intersection(wodu_range_bli, northamerica) %>% summarise()
rndu_range_bli_NA <- st_intersection(rndu_range_bli, northamerica) %>% summarise()

# Clip breeding ranges down to the Miss and Central Flyway
mall_bli_assignment_range <- st_intersection(mall_range_bli_NA, ms_cn_flyway) %>%
  summarise()
mapview(mall_bli_assignment_range)

wodu_bli_assignment_range <- st_intersection(wodu_range_bli_NA, ms_cn_flyway) %>%
  summarise()
mapview(wodu_bli_assignment_range)

rndu_bli_assignment_range <- st_intersection(rndu_range_bli_NA, ms_cn_flyway) %>%
  summarise()
mapview(rndu_bli_assignment_range)

# save assignment ranges
write_rds(mall_bli_assignment_range, "data/processed/spp_assignment_ranges/mall_bli_assignment_range.rds")
write_rds(wodu_bli_assignment_range, "data/processed/spp_assignment_ranges/wodu_bli_assignment_range.rds")
write_rds(rndu_bli_assignment_range, "data/processed/spp_assignment_ranges/rndu_bli_assignment_range.rds")

############ 
## Figures to compare Bird life international ranges to ebird ranges

# good mapping projection
crs <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# states and provinces
northamerica <- ne_states(iso_a2 = c("CA", "US"), returnclass = "sf") %>% 
  filter(!name %in% c("Alaska", "Hawaii")) %>% 
  st_transform(crs = crs)

# transform objects for mapping
ms_cn_flyway_proj <- ms_cn_flyway %>% 
  st_transform(crs = crs)

mall_ebird_assignment_range_proj <- mall_ebird_assignment_range %>% 
  st_transform(crs = crs)

mall_bli_assignment_range_proj <- mall_bli_assignment_range %>% 
  st_transform(crs = crs)

wodu_ebird_assignment_range_proj <- wodu_ebird_assignment_range %>% 
  st_transform(crs = crs)

wodu_bli_assignment_range_proj <- wodu_bli_assignment_range %>% 
  st_transform(crs = crs)

rndu_ebird_assignment_range_proj <- rndu_ebird_assignment_range %>% 
  st_transform(crs = crs)

rndu_bli_assignment_range_proj <- rndu_bli_assignment_range %>% 
  st_transform(crs = crs)

# mall figure

mall_range_comp <- ggplot() +
  geom_sf(data = ms_cn_flyway_proj, fill = "lightgray") +
  geom_sf(data = northamerica, fill = NA) +
  geom_sf(data = mall_bli_assignment_range_proj, fill = "red", alpha = 0.4, color = "black") + 
  geom_sf(data = mall_ebird_assignment_range_proj, fill = "blue", alpha = 0.4, color = "black") +
  theme_minimal() +
  labs(title = "Mallard range source comparisons",
       subtitle = "Pink = BirdLife International | Blue = eBird | Purple = overlap\n(MS and CN flyways in gray)")
mall_range_comp
ggsave(here::here('figures/cal_equation_selection/mall_range_comparison.png'), height = 6, width = 8, units = 'in', dpi = 500)

# wodu

wodu_range_comp <- ggplot() +
  geom_sf(data = ms_cn_flyway_proj, fill = "lightgray") +
  geom_sf(data = northamerica, fill = NA) +
  geom_sf(data = wodu_bli_assignment_range_proj, fill = "red", alpha = 0.4, color = "black") + 
  geom_sf(data = wodu_ebird_assignment_range_proj, fill = "blue", alpha = 0.4, color = "black") +
  theme_minimal() +
  labs(title = "Wood duck range source comparisons",
       subtitle = "Pink = BirdLife International | Blue = eBird | Purple = overlap\n(MS and CN flyways in gray)")
wodu_range_comp
ggsave(here::here('figures/cal_equation_selection/wodu_range_comparison.png'), height = 6, width = 8, units = 'in', dpi = 500)

# rndu figure

rndu_range_comp <- ggplot() +
  geom_sf(data = ms_cn_flyway_proj, fill = "lightgray") +
  geom_sf(data = northamerica, fill = NA) +
  geom_sf(data = rndu_bli_assignment_range_proj, fill = "red", alpha = 0.4, color = "black") + 
  geom_sf(data = rndu_ebird_assignment_range_proj, fill = "blue", alpha = 0.4, color = "black") +
  theme_minimal() +
  labs(title = "Ring-necked duck range source comparisons",
       subtitle = "Pink = BirdLife International | Blue = eBird | Purple = overlap\n(MS and CN flyways in gray)")
rndu_range_comp
ggsave(here::here('figures/cal_equation_selection/rndu_range_comparison.png'), height = 6, width = 8, units = 'in', dpi = 500)

