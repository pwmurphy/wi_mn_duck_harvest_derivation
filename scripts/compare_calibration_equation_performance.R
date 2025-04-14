### compare GS and MA calibration equations (Kusack et al. 2023) performance on all local reference samples #####

library(isocat)
library(raster)
library(terra)
library(sf)
library(tidyverse)
library(mapview)
library(readr)
library(assignR)
library(RColorBrewer)
library(rnaturalearth)

## prepare environment
rm(list=ls()) # clear the global environment to save on space
# rm(object_name) # to remove specific objects
set.seed(725)
dev.off()
select <- dplyr::select
custom.palette <- colorRampPalette(brewer.pal(9, "YlGnBu")) # Color palette
gc()

## Load in data needed to run the spatial assignment
# dataframe with sample d2h values and metadata
spp_d2h <- read_rds("data/processed/cleaned_feather_results/reference_data_df.rds") %>% 
  filter(age == "local") %>% # for now, not using AHY or HY birds 
  filter(!(state == "WI" & species == "RNDU")) %>% # the two rndu samples from WI fall outside eBird breeding range crop
  mutate(species = tolower(species))

# dataframe with calibration equation values from Jackson's paper
cal_equations <- tibble(
  precip_iso = c("gs", "ma", "gs", "ma"),
  forage_guild = c("dabbler", "dabbler", "diver", "diver"),
  slope = c(0.7, 0.6, 0.5, 0.5),
  intercept = c(-69.9, -71.7, -82.6, -78.4),
  eq_sd = c(17.7, 17.2, 14.7, 14.1)
)

# species eBird breeding range maps to crop the assignment range to
files <- list.files("data/processed/spp_assignment_ranges/", pattern = "ebird", full.names = TRUE) # find range map files

spp_assignment_ranges <- map_dfr(files, ~ { # read in all files and combine to one object
  species <- substr(basename(.x), 1, 4)  
  read_rds(.x) %>% mutate(species = species, .before = geom)
}) %>%
  mutate(forage_guild = case_when( # specify species foraging guild (so we know which calibration equation to use)
    species %in% c("wodu", "mall") ~ "dabbler",
    species == "rndu" ~ "diver"
  ), .after = species) %>% 
  mutate(common_name = case_when( # add in a column for common names (to pull later for use in figure titles)
    species == "mall" ~ "Mallard",
    species == "wodu" ~ "Wood Duck",
    species == "rndu" ~ "Ring-necked Duck"
  ), .before = geom)

## Load in precipitation isoscapes and associated error isoscapes

# growing-season isoscape
isoscape_gs <- rast("data/raw/precipitation_isoscapes/gs_isoscape/d2h_GS.tif")
crs(isoscape_gs) <- "EPSG:4326"

isoscape_gs_se <- rast("data/raw/precipitation_isoscapes/gs_isoscape/d2h_se_GS.tif")
crs(isoscape_gs_se) <- "EPSG:4326"

# mean annual isoscape
isoscape_ma <- rast("data/raw/precipitation_isoscapes/ma_isoscape/d2h_MA.tif")
crs(isoscape_ma) <- "EPSG:4326"

isoscape_ma_se <- rast("data/raw/precipitation_isoscapes/ma_isoscape/d2h_se_MA.tif")
crs(isoscape_ma_se) <- "EPSG:4326"

## load up other data for mapping and specify mapping preferences
# preferred mapping crs
crs <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # set nice projection for mapping

# flyways, countries, state/province outlines
northamerica <- ne_states(iso_a2 = c("CA", "US"), returnclass = "sf") %>% 
  filter(!name %in% c("Alaska", "Hawaii")) %>% 
  st_transform(crs = crs)

ms_cn_flyway_proj <- read_rds("data/processed/isolated_flyways/ms_cn_contin_flyway.rds") %>% 
  st_transform(crs = crs)

wi <- ne_states(country = 'united states of america', returnclass = 'sf') %>%
  filter(postal == 'WI')

wi <- wi %>%
  st_transform(., st_crs(4326)) %>%
  dplyr::select(postal)

mn <- ne_states(country = 'united states of america', returnclass = 'sf') %>%
  filter(postal == 'MN')

mn <- mn %>%
  st_transform(., st_crs(4326)) %>%
  dplyr::select(postal)

wi_mn <- rbind(wi, mn) # combine Wi and MN so we can add them to maps as one thing

#################################################

## perform the spatial assignment for all reference samples and check the accuracy of the assignment 

# below: function to run through samples from one species at a time and do spatial assignments with both MA and GS calibration equations
# (1) calibrate and crop isoscapes for GS and MA
# (2) perform spatial assignment with isocat function
# (3) get binary 2:1 odds assignments with AssignR
# (4) make nice ggplots of the species combined 2:1 odds origins
# (5) calculate the number that were actually correctly assigned by MA and GS equations
# (6) make nice ggplots comparing how well the equations did on predicting reference sample origins by mapping out actual origin locations

check_origins <- function(sp_to_run) {
  
  # filter spp_assignment_ranges df
  this_species <- spp_assignment_ranges %>% filter(species == sp_to_run)
  
  # Extract species-specific info from the row
  species <- sp_to_run
  foraging_guild <- this_species$forage_guild
  common_name <- this_species$common_name
  
  # Calibrate and crop isoscapes for growing season (gs) and mean annual (ma)
  # gs
  cal_eq_gs <- cal_equations %>% filter(foraging_guild == forage_guild & precip_iso == "gs")
  
  cal_iso_gs <- (isoscape_gs * cal_eq_gs$slope + cal_eq_gs$intercept) %>% 
    mask(this_species) %>%  
    crop(ext(this_species))
  
  cal_iso_gs_se <- (isoscape_gs_se * cal_eq_gs$slope) %>% 
    mask(this_species) %>%  
    crop(ext(this_species))
  
  # ma 
  cal_eq_ma <- cal_equations %>% filter(foraging_guild == forage_guild & precip_iso == "ma")
  
  cal_iso_ma <- (isoscape_ma * cal_eq_ma$slope + cal_eq_ma$intercept) %>% 
    mask(this_species) %>%  
    crop(ext(this_species))
  
  cal_iso_ma_se <- (isoscape_ma_se * cal_eq_ma$slope) %>% 
    mask(this_species) %>%  
    crop(ext(this_species))
  
  ## Spatial assignment w isocat
  sp_d2h <- spp_d2h %>% filter(species == sp_to_run)
  
  origins_gs <- isotopeAssignmentModel(ID = sp_d2h$cornell_sample_id, 
                                       isotopeValue = sp_d2h$d2h_vs_vsmow_original, 
                                       SD_indv = cal_eq_gs$eq_sd, 
                                       precip_raster = raster(cal_iso_gs), 
                                       precip_SD_raster = raster(cal_iso_gs_se))  # Output is a large rasterstack
  
  origins_ma <- isotopeAssignmentModel(ID = sp_d2h$cornell_sample_id, 
                                       isotopeValue = sp_d2h$d2h_vs_vsmow_original, 
                                       SD_indv = cal_eq_ma$eq_sd, 
                                       precip_raster = raster(cal_iso_ma), 
                                       precip_SD_raster = raster(cal_iso_ma_se))  # Output is a large rasterstack
  
  ## Binary assignment w AssignR
  bin_origins_gs <- assignR::qtlRaster(rast(origins_gs), threshold = 0.66, thresholdType = "prob") * 1 
  names(bin_origins_gs) <- rep("d2h_gs", length(names(bin_origins_gs)))
  bin_all_origins_gs <- app(bin_origins_gs, sum)  # get combined likely origins map across all birds
  
  bin_origins_ma <- assignR::qtlRaster(rast(origins_ma), threshold = 0.66, thresholdType = "prob") * 1
  names(bin_origins_ma) <- rep("d2h_ma", length(names(bin_origins_ma)))
  bin_all_origins_ma <- app(bin_origins_ma, sum)  # get combined likely origins map across all birds
  
  # prepare species origin rasters for plotting
  bin_all_origins_gs_df <- bin_all_origins_gs %>% 
    project(crs, method = 'near') %>% 
    as.data.frame(xy = TRUE)
  
  bin_all_origins_ma_df <- bin_all_origins_ma %>% 
    project(crs, method = 'near') %>% 
    as.data.frame(xy = TRUE)
  
  # plot combined 2:1 likely origins maps
  gs_combined_origins <- ggplot() +
    geom_sf(data = ms_cn_flyway_proj, fill = "lightgray") +
    geom_tile(data = bin_all_origins_gs_df, aes(x = x, y = y, fill = sum)) +
    geom_sf(data = northamerica, fill = NA) +
    scale_fill_gradientn(colors = rev(custom.palette(150)), name = " ") +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(title = paste0("Likely origins of all local ", this_species$common_name, " reference samples\nunder 2:1 odds ratio, n = ", nrow(sp_d2h), "\n(using growing season calibration)"), x = " ", y = " ")
  
  print(gs_combined_origins)
  gs_file_path <- here::here(paste0('figures/cal_equation_selection/local_', species, '_ref_origins_gs.png'))
  ggsave(gs_file_path, gs_combined_origins, height = 6, width = 8, units = 'in', dpi = 500)
  message("Saved growing season combined origins figure for ", species, " to: ", gs_file_path)
  
  ma_combined_origins <- ggplot() +
    geom_sf(data = ms_cn_flyway_proj, fill = "lightgray") +
    geom_tile(data = bin_all_origins_ma_df, aes(x = x, y = y, fill = sum)) +
    geom_sf(data = northamerica, fill = NA) +
    scale_fill_gradientn(colors = rev(custom.palette(150)), name = " ") +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(title = paste0("Likely origins of all local ", this_species$common_name, " reference samples\nunder 2:1 odds ratio, n = ", nrow(sp_d2h), "\n(using mean annual calibration)"), x = " ", y = " ")
  
  print(ma_combined_origins)
  ma_file_path <- here::here(paste0('figures/cal_equation_selection/local_', species, '_ref_origins_ma.png'))
  ggsave(ma_file_path, ma_combined_origins, height = 6, width = 8, units = 'in', dpi = 500)
  message("Saved mean annual combined origins figure for ", species, " to: ", gs_file_path)
  
  ## Check correctness of assignments
  birds <- as.vector(sp_d2h$d2h_vs_vsmow_original)  # all local birds of this species
  sp_sf <- sp_d2h %>%
    st_as_sf(coords = c('longitude', 'latitude'), crs = 4326)
  
  # Convert sf points to a SpatVector
  sp_vect <- vect(sp_sf)
  
  # run over all birds to get raster values
  correct_origins <- map_dfr(seq_along(birds), function(j) {
    raster_value_gs <- terra::extract(bin_origins_gs[[j]], sp_vect[j, ])
    raster_value_ma <- terra::extract(bin_origins_ma[[j]], sp_vect[j, ])
    
    tibble(
      id = sp_sf$cornell_sample_id[j],
      species = species,
      new_sample = sp_d2h$new_sample_2025[j],
      latitude = sp_d2h$latitude[j],
      longitude = sp_d2h$longitude[j],
      gs_correct = raster_value_gs[, 2],
      ma_correct = raster_value_ma[, 2], 
      correct = case_when(
        gs_correct == 1 & ma_correct == 1 ~ "Correct by both eq.",
        gs_correct == 1 & ma_correct == 0 ~ "GS only correct",
        gs_correct == 0 & ma_correct == 1 ~ "MA only correct",
        gs_correct == 0 & ma_correct == 0 ~ "Neither eq. correct"))
  })
  
  # map out the origin locations and whether they were correctly identified by the two equations
    custom_colors <- c(
      "Correct by both eq." = "#2fb52c",
      "GS only correct" = "#eca721",
      "MA only correct" = "#d8c53b",
      "Neither eq. correct" = "#3e50f7"
    )
    
  correct_plot <- ggplot() +
    geom_jitter(data = correct_origins, aes(x = longitude, y = latitude, color = correct), 
                width = 0.2, height = 0.2, alpha = 0.7, size = 2.5) +  # jitter points
    geom_sf(data = wi_mn, fill = NA) +  # background map
    scale_color_manual(values = custom_colors) + 
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs( y = "", x = "", color = " ") + 
    ggtitle(paste0("Comparing Calibration Equation Performance: \n All Local ", common_name, " Reference Samples, n = ", nrow(correct_origins)))
  
  print(correct_plot)
  correct_file_path <- here::here(paste0('figures/cal_equation_selection/comp_local_', species, '_origins.png'))
  ggsave(correct_file_path, correct_plot, height = 6, width = 8, units = 'in', dpi = 500)
  message("Saved equation comparison figure for ", species, " to: ", correct_file_path)
  
  # filter to only show the newly processed samples from 2024/2025 batches
  correct_origins_new <-correct_origins %>% filter(new_sample == 1)
  
  new_correct_plot <- ggplot() +
    geom_jitter(data = correct_origins_new, aes(x = longitude, y = latitude, color = correct), 
                width = 0.2, height = 0.2, alpha = 0.7, size = 2.5) +  # jitter points
    geom_sf(data = wi_mn, fill = NA) +  # background map
    scale_color_manual(values = custom_colors) + 
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs( y = "", x = "", color = " ") + 
    ggtitle(paste0("Comparing Calibration Equation Performance: \n New Local ", common_name, " Reference Samples, n = ", nrow(correct_origins_new)))
  
  print(new_correct_plot)
  new_correct_file_path <- here::here(paste0('figures/cal_equation_selection/comp_new_local_', species, '_origins.png'))
  ggsave(new_correct_file_path, new_correct_plot, height = 6, width = 8, units = 'in', dpi = 500)
  message("Saved new sample equation comparison figure for ", species, " to: ", new_correct_file_path)
  
  # Save the output to the environment
  assign(paste0("correct_origins_", species), correct_origins, envir = .GlobalEnv)
  assign(paste0("origins_gs_", species), origins_gs, envir = .GlobalEnv)
  assign(paste0("origins_ma_", species), origins_ma, envir = .GlobalEnv)
  assign(paste0("bin_origins_gs_", species), bin_origins_gs, envir = .GlobalEnv)
  assign(paste0("bin_origins_ma_", species), bin_origins_ma, envir = .GlobalEnv)
  assign(paste0("bin_all_origins_gs_", species), bin_all_origins_gs, envir = .GlobalEnv)
  assign(paste0("bin_all_origins_ma_", species), bin_all_origins_ma, envir = .GlobalEnv)
}

# run for each species, look at output (may take a few minutes when there are many samples)
check_origins("wodu")
check_origins("rndu")
check_origins("mall")

# combine origin accuracy output into one table 
cal_eq_performance_tbl <- rbind(correct_origins_mall, correct_origins_wodu, correct_origins_rndu) 

# save table of results for each sample
write_csv(cal_eq_performance_tbl, "output/cal_equation_selection/eq_performance_by_sample.csv")

## summarize origin accuracy results by species
origins_summary_all <- cal_eq_performance_tbl %>% 
  select(!correct) %>% 
  pivot_longer(c(gs_correct, ma_correct), names_to = "equation", values_to = "correct") %>% 
  group_by(species, equation) %>%
  summarise(
    total_correct = sum(correct, na.rm = TRUE),
    total_all = n(),  # Total number of samples per group
    .groups = "drop"
  ) %>%
  mutate(
    total_percentage = round((total_correct / total_all) * 100),
    id = paste(species, equation, sep = "_")
  ) %>%
  select(-species, -equation)

# also, separately reporting the subset of values from newly processed samples because these were not included in the datasets used to create the calibration equation (out of sample performance indicators)
origins_summary_tbl <- cal_eq_performance_tbl %>% 
  select(!correct) %>% 
  pivot_longer(c(gs_correct, ma_correct), names_to = "equation", values_to = "correct") %>% 
  group_by(species, equation, new_sample) %>%
  summarise(
    new_correct = sum(correct, na.rm = TRUE),
    new_all = n(),  # Total number of new samples per group
    .groups = "drop"
  ) %>%
  mutate(
    new_percentage = round((new_correct / new_all) * 100),
    id = paste(species, equation, sep = "_")
  ) %>%
  filter(new_sample == 1) %>%
  right_join(origins_summary_all, by = "id") %>%
  mutate(
    new_samples = paste0(new_correct, "/", new_all, " (", new_percentage, "%)"),
    all_samples = paste0(total_correct, "/", total_all, " (", total_percentage, "%)")
  ) %>%
  select(species, equation, all_samples, new_samples) %>% 
  mutate(order = c(1,2,5,6,3,4)) %>% 
  arrange(order) %>% 
  select(!order)

# save the output 
write_csv(origins_summary_tbl, "output/cal_equation_selection/correct_origins_local_summary_tbl.csv")

######################################

# example plots of two bird origins for 3/25 meeting-- one correctly assigned and one not

## extract origin plots from stack
# probability surface
WI_23_83 <- origins_gs_wodu[[26]]
WI_23_9 <- origins_gs_wodu[[21]]

# binary surface
WI_23_83_bin <- bin_origins_gs_wodu[[26]]
WI_23_9_bin <- bin_origins_gs_wodu[[21]]

# Convert rasters to df for ggplot
WI_23_83_df <- WI_23_83 %>% 
  rast() %>% 
  project(crs, method = 'near') %>% 
  as.data.frame(xy = TRUE)

WI_23_9_df <- WI_23_9 %>% 
  rast() %>% 
  project(crs, method = 'near') %>% 
  as.data.frame(xy = TRUE)

WI_23_83_bin_df <- WI_23_83_bin %>% 
  project(crs, method = 'near') %>% 
  as.data.frame(xy = TRUE)

WI_23_9_bin_df <- WI_23_9_bin %>% 
  project(crs, method = 'near') %>% 
  as.data.frame(xy = TRUE)

# get an sf object of the actual origin location of the reference sample
WI_23_83_pt <- spp_d2h %>% 
  filter(cornell_sample_id == "WI-23-83") %>% 
  select(cornell_sample_id, latitude, longitude) %>% 
  st_as_sf(coords = c('longitude', 'latitude'), crs = 4326) %>% 
  st_transform(crs = crs)

WI_23_9_pt <- spp_d2h %>% 
  filter(cornell_sample_id == "WI-23-9") %>% 
  select(cornell_sample_id, latitude, longitude) %>% 
  st_as_sf(coords = c('longitude', 'latitude'), crs = 4326) %>% 
  st_transform(crs = crs)

## plot the full probability surface for each sample
# plot 23-83 (correctly identified origin)
WI_23_83_plot <- ggplot() +
  geom_sf(data = ms_cn_flyway_proj, fill = "lightgray") +
  geom_tile(data = WI_23_83_df, aes(x = x, y = y, fill = WI.23.83)) +  # Convert to factor
  geom_sf(data = northamerica, fill = NA) +
  geom_sf(data = WI_23_83_pt, color = "red") +
  scale_fill_gradientn(colors = rev(custom.palette(150)), name = "Origin \nProbability") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Origin probability of wood duck reference sample WI-23-83 \n(using growing season precip. calibration)", 
       x = " ", y = " ")

WI_23_83_plot
ggsave(here::here('figures/cal_equation_selection/example_wodu_WI2383.png'), height = 6, width = 8, units = 'in', dpi = 500)

# plot 23-9 (incorrectly identified origin)
WI_23_9_plot <- ggplot() +
  geom_sf(data = ms_cn_flyway_proj, fill = "lightgray") +
  geom_tile(data = WI_23_9_df, aes(x = x, y = y, fill = WI.23.9)) +  # Convert to factor
  geom_sf(data = northamerica, fill = NA) +
  geom_sf(data = WI_23_9_pt, color = "red") +
  scale_fill_gradientn(colors = rev(custom.palette(150)), name = "Origin \nProbability") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Origin probability of wood duck reference sample WI-23-9 \n(using growing season precip. calibration)", 
       x = " ", y = " ")

WI_23_9_plot
ggsave(here::here('figures/cal_equation_selection/example_wodu_WI239.png'), height = 6, width = 8, units = 'in', dpi = 500)

## plot the binary surface for each sample
# plot 23-83 (correctly identified origin)
WI_23_83_bin_plot <- ggplot() +
  geom_sf(data = ms_cn_flyway_proj, fill = "lightgray") +
  geom_tile(data = WI_23_83_bin_df, aes(x = x, y = y, fill = as.factor(d2h_gs))) +  # Convert to factor
  geom_sf(data = northamerica, fill = NA) +
  geom_sf(data = WI_23_83_pt, color = "red") +
  scale_fill_manual(
    values = c("0" = "navy", "1" = "lightblue"),  # Assign colors
    name = "Origin Likelihood",  # Legend title
    labels = c("0" = "Not likely origin", "1" = "66% likely origin")  # Rename categories
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Likely origins of wood duck reference sample WI-23-83 \nunder 2:1 odds ratio \n(using growing season precip. calibration)", 
       x = " ", y = " ")

WI_23_83_bin_plot
ggsave(here::here('figures/cal_equation_selection/example_wodu_binary_WI2383.png'), height = 6, width = 8, units = 'in', dpi = 500)

# plot 23-9 (incorrectly identified origin)
WI_23_9_bin_plot <- ggplot() +
  geom_sf(data = ms_cn_flyway_proj, fill = "lightgray") +
  geom_tile(data = WI_23_9_bin_df, aes(x = x, y = y, fill = as.factor(d2h_gs))) +  # Convert to factor
  geom_sf(data = northamerica, fill = NA) +
  geom_sf(data = WI_23_9_pt, color = "red") +
  scale_fill_manual(
    values = c("0" = "navy", "1" = "lightblue"),  # Assign colors
    name = "Origin Likelihood",  # Legend title
    labels = c("0" = "Not likely origin", "1" = "66% likely origin")  # Rename categories
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Likely origins of wood duck reference sample WI-23-9 \nunder 2:1 odds ratio \n(using growing season precip. calibration)", 
       x = " ", y = " ")

WI_23_9_bin_plot
ggsave(here::here('figures/cal_equation_selection/example_wodu_binary_WI239.png'), height = 6, width = 8, units = 'in', dpi = 500)

# 
# ################ looking at other summary functions in isocat
# # cumulative sum surface 
# CumSum <- makecumsumSurface(prob[[10]])
# #(Each value represents the cumulative sum of all values ≤ that cell.)
# plot(CumSum)
# 
# Odds <- makeOddsSurfaces(prob[[10]])
# # Converts probabilities into odds
# plot(Odds)
