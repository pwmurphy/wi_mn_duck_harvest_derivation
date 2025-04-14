######### preparing isotope results by combine with collection data #######
######## for all species, all sample types, from MN DNR and WI DNR ########

# libraries
library(fs)
library(readxl)
library(tidyverse)
library(janitor)
library(arcpullr)
library(tigris) # for mn county data
options(tigris_use_cache = TRUE)

# load spatial data (need to join counties to WI reference feathers -----------------
# load county shapefile (using arcpullr from WDNR server)
wdnr_server <-"https://dnrmaps.wi.gov/arcgis/rest/services/"
counties <- "DW_Map_Dynamic/EN_Basic_Basemap_WTM_Ext_Dynamic_L16/MapServer/3"
wi_counties_url <- paste(wdnr_server, counties, sep ="/")
wi_counties <- get_spatial_layer(wi_counties_url) %>% select(county = COUNTY_NAME)


# locate sample IDs for submitted feathers --------------------------------

# find Cornell submission sheets
feathers_submitted_df <- dir_ls(
  here::here('data/raw/feather_collection_data/submission_sheets') 
) %>%
  map_dfr(., read_excel) %>% # read, combine excel worksheets
  select(cornell_sample_id = WI_ID_Code)
feathers_submitted_df


# now join samples with feather collection metadata -----------------------

# WINGBEE

# read in wingbee samples (2019 and 2020)
wingbee_feathers_2019_2020_df <- dir_ls(
  here::here('data/raw/feather_collection_data'), 
  glob = '*wingbee_feathers.xlsx'
) %>%
  map_dfr(., read_excel) # read and combine
wingbee_feathers_2019_2020_df

# clean up 2019, 2020 wingbee feathers
wingbee_feathers_2019_2020_df <- wingbee_feathers_2019_2020_df %>%
  mutate(
    harvest_date = ymd(harvest_date), # fix dates
    across(where(is.character), ~na_if(., ".")), # replace . with NA
    source = 'wingbee', # add source column
    # fix species ID
    species = case_when(
      species == 'Mall' ~ 'MALL',
      TRUE ~ species
    )
  ) %>%
  rename(
    collection_date = harvest_date,
    cornell_sample_id = wi_sample_id
  )
wingbee_feathers_2019_2020_df

# now read in feathers from 2021 wingbee (samples submitted in March 2022)
# this is the full list of data we received from stephen chandler;
# we only sent in a subset of these to cornell
wingbee_feathers_2021_df <- read_csv(here::here('data/raw/feather_collection_data/complete_list_2021_wingbee_feathers.csv'))

# clean up
wingbee_feathers_2021_df <- wingbee_feathers_2021_df %>%
  mutate(
    cornell_sample_id = fws_id, # in 2022, we just used this rather than making new IDs
    source = 'wingbee', # add source column
    state = 'WI',
    feather = 'P1',
    town = as.character(NA),
    comments = as.character(NA),
    collection_date = mdy(collection_date)
  ) %>%
  select(
    cornell_sample_id, 
    fws_id, 
    species, 
    state, 
    county, 
    town, 
    collection_date, 
    age, 
    sex, 
    feather, 
    comments, 
    source
  )
wingbee_feathers_2021_df

# combine all wingbee metadata
wingbee_feathers_df <- wingbee_feathers_2019_2020_df %>%
  bind_rows(., wingbee_feathers_2021_df)
wingbee_feathers_df

# this contains extra data (feathers not included in March 2022 submission)
wingbee_feathers_df %>% print(n=Inf)

# REFERENCE

# manually looked up Lat Longs for the 2023 feathers because they weren't listed 
# generally made the coordinate the middle of the waterbody unless more specific info given or had coordinates from previous years

# now read in reference samples
reference_feathers_df <- dir_ls(
  here::here('data/raw/feather_collection_data'),
  glob = '*reference_feathers.xlsx'
) %>%
  map_dfr(., read_excel)
reference_feathers_df

# more cleaning
reference_feathers_df <- reference_feathers_df %>%
  select(-easting, -northing, -...17) %>% # drop these
  drop_na(latitude) %>% # remove entries with no lat/long 
  mutate(
    collection_date = ymd(collection_date), # fix dates
    across(where(is.character), ~na_if(., ".")), # replace . with NA
    source = 'reference'
  ) %>%
  rename(cornell_sample_id = wi_sample_id)
reference_feathers_df

# now convert to sf object to extract county info
reference_feathers_sf <- reference_feathers_df %>%
  st_as_sf(coords = c('longitude', 'latitude'), crs = 4326)
reference_feathers_sf

# extract county for each point
reference_feathers_sf <- st_join(reference_feathers_sf, wi_counties, join = st_within) %>%
  select(cornell_sample_id, county) %>%
  st_drop_geometry()
reference_feathers_sf

# now join this back to reference feathers df to add counties
reference_feathers_df <- reference_feathers_df %>%
  left_join(., reference_feathers_sf)
reference_feathers_df

# HUNTER (2021)
hunter_feathers_2021_df <- read_excel(here::here('data/raw/feather_collection_data/rndu_hunter_harvested_wings_2022.xlsx'))

hunter_feathers_2021_df <- hunter_feathers_2021_df %>%
  rename(cornell_sample_id = wi_sample_id) %>%
  mutate(
    collection_date = ymd(collection_date),
    county = str_to_title(county),
    feather = 'P1',
    source = 'hunter',
    state = 'WI'
  )
hunter_feathers_2021_df

# PUT ALL OF IT TOGETHER

# combine with wingbee data
all_feathers_metadata_df <- wingbee_feathers_df %>% # wingbee
  bind_rows(., reference_feathers_df) %>% # reference
  bind_rows(., hunter_feathers_2021_df) %>% # hunter
  # only keep sample IDs here that match the ones that were submitted
  filter(cornell_sample_id %in% c(feathers_submitted_df %>% pull()))
all_feathers_metadata_df


# clean isotope results from cornell -------------------------------------------------------

# batches to loop through and clean up
batches <- c(
  'batch1', 
  'batch2', 
  'batch3',
  'batch4'
)

# create empty list to store cleaned data frames
cornell_results_df <- list()

# loop through each batch
for (i in seq_along(batches)) {
  
  # batch i excel file location
  file <- str_glue("data/raw/cornell_coil_data/cornell_d2H_analysis_results_{batches[i]}.xlsx")
  
  # find sheet names in this file
  sheets <- excel_sheets(file)
  sheets <- sheets %>%
    str_subset(pattern = "run") # discard cover letter sheet
  
  # combine all sheets (run 1, 2, ...)
  df <- map_dfr(sheets, ~read_excel(file, sheet = .x, skip = 8), .id = "run_id") %>%
    clean_names()
  
  # fix a bit
  df <- df %>%
    rename(
      d2h_vs_vsmow_original = d2h_vs_vsmow_5,
      d2h_vs_vsmow_2017 = d2h_vs_vsmow_6,
      cornell_sample_id = sample_id
    ) %>%
    # add run name
    mutate(run_id = str_c('run_', run_id))
  df
  
  # find all cornell standards across runs
  cornell_standards_df <- df %>%
    filter(str_detect(cornell_sample_id, "std")) %>%
    mutate(
      across(weight_mg:d2h_vs_vsmow_2017, as.numeric),
      sample_type = 'standard'
    )
  cornell_standards_df
  
  # find all cornell samples across runs
  # we don't want these 'samples'
  non_samples <- "std|Quality Control Data|Benzoic Acid|In-house standards used for precision purposes|Sample ID|CBS|KHS|Keratin|CBC|In-house standards used for percent element calculation"
  
  cornell_samples_df <- df %>%
    filter(!str_detect(cornell_sample_id, non_samples)) %>%
    mutate(across(weight_mg:d2h_vs_vsmow_2017, as.numeric)) %>%
    # remove spaces, fix one bad name
    mutate(
      cornell_sample_id = str_remove(cornell_sample_id, ' '),
      cornell_sample_id = case_when(
        # from second batch
        cornell_sample_id == 'LM2004' ~ 'LM20004',
        # from third batch
        cornell_sample_id == 'M5M2694' ~ 'M552694',
        cornell_sample_id == 'M519194' ~ 'M516194',
        cornell_sample_id == 'M640640' ~ 'M610640',
        cornell_sample_id == 'MM52695' ~ 'M552695',
        cornell_sample_id == 'M519196' ~ 'M516196',
        # otherwise they are ok as is
        TRUE ~ cornell_sample_id
      ),
      sample_type = 'feather'
    )
  cornell_samples_df
  
  # join all cornell data back together
  cornell_df <- cornell_samples_df %>%
    bind_rows(., cornell_standards_df) %>%
    select(
      run_id, 
      cornell_sample_id, 
      weight_mg, 
      h2_amp, 
      percent_h, 
      d2h_vs_vsmow_original, 
      d2h_vs_vsmow_2017,
      sample_type
    ) %>%
    # keep batch id
    mutate(batch_id = str_glue({batches[i]}))
  cornell_df
  
  # store data frame in list
  cornell_results_df[[i]] <- cornell_df
  
}

# combine each batch
cornell_results_df <- cornell_results_df %>%
  bind_rows() # list of data frames to single data frame
cornell_results_df

# link back to sample metadata --------------------------------------------

feather_results_df <- all_feathers_metadata_df %>%
  left_join(., cornell_results_df) %>%
  select(cornell_sample_id, fws_id, batch_id, run_id, source, sample_type, collection_date, species, age, sex, feather, state, county, town, latitude, longitude, weight_mg:d2h_vs_vsmow_2017, comments) %>% 
  filter(!cornell_sample_id %in% c("LR21020", "M621633", "M598489"))
feather_results_df 
# missing results for LR21020 -- not in Cornell sheets. Also two Wingbee feathers did not get d2h results (low sample wt or other issue): M621633 and M598489

feather_results_df %>%
  filter(sample_type == 'feather') %>%
  print(n=Inf)

# fix incorrect dates for two entries
feather_results_df$collection_date[feather_results_df$cornell_sample_id == 'R20071'] <- "2020-10-31"
feather_results_df$collection_date[feather_results_df$cornell_sample_id == 'LW21002'] <- "2021-07-14"

# create 'year' column 
feather_results_df <- feather_results_df %>% 
  mutate(year = year(collection_date))
feather_results_df

unique(feather_results_df$sex) # one sample from 2023 with sex mislabeled

# 8-4-22 Read in MN reference feathers provided by Bruce Davis 
mn_reference_feathers_2020 <- read_csv(
  "data/raw/mn_dnr_data/mn_reference_feathers_2020.csv",
  col_types = cols(
    fws_id = col_character(),
    batch_id = col_character(),
    feather = col_character(),
    latitude = col_number(),
    longitude = col_number(), 
    d2h_vs_vsmow_2017 = col_number()
  )
) %>%
  mutate(collection_date = mdy(collection_date))
mn_reference_feathers_2020
table(mn_reference_feathers_2020$species)

# clean up mn_reference_feathers2020
mn_reference_feathers_2020 <- mn_reference_feathers_2020 %>%
  mutate(
    # fix species ID
    age = case_when(
      age == 'HY' ~ 'hy',
      age == 'AHY' ~ 'ahy',
      age == 'L' ~ 'local',
      TRUE ~ age
    ),
    sex = case_when(
      sex == 'M' ~ 'male',
      sex == 'F' ~ 'female',
      TRUE ~ sex
    ),
    sample_type = case_when(
      is.na(sample_type) ~ 'feather',
      TRUE ~ sample_type
    ),
  ) %>%
  mutate(
    run_id = case_when(
      run_id == '1' ~ 'run_1',
      run_id == '2' ~ 'run_2')
  )
mn_reference_feathers_2020

# merge MN reference feather to all feathers data frame
feather_results_df <- rbind(feather_results_df, mn_reference_feathers_2020) %>%
  mutate(doy = yday(collection_date)) 
feather_results_df

###NOTE - 6-22-22:  Need to change samples that have "volunteer submissions" in "comments' row to "hunter" under "source" column.  Right now, currently labeled as "wingbee"
unique(feather_results_df$comments)

# check these
feather_results_df %>%
  filter(comments == 'voluneteer submission') %>%
  print(n=Inf)

# fix
feather_results_df <- feather_results_df %>%
  # change source to hunter if it is a volunteer submission;
  # otherwise keep source the same as listed
  mutate(
    source = case_when(
      comments == 'voluneteer submission' ~ 'hunter',
      TRUE ~ source
    )
  )
feather_results_df

# fix state abbreviation (there's a Wi that should be WI)
feather_results_df <- feather_results_df %>%
  mutate(
    state = case_when(
      state == 'Wi' ~ 'WI',
      TRUE ~ state
    )
  )
feather_results_df

# clean up, keep relevant stuff
feather_results_df <- feather_results_df %>%
  select(
    species, 
    cornell_sample_id, 
    batch_id, 
    run_id, 
    source, 
    sample_type, 
    collection_date, 
    year, 
    doy, 
    age, 
    sex, 
    feather, 
    state, 
    county,
    town,
    latitude, 
    longitude, 
    weight_mg:comments
    )

feather_results_df

## 

# 3-28-25 add in all newly received MN data-- reference and wingbee
# Bruce sent it to us w same format as our feather_results_df so minimal cleaning needed

mn_feathers_2025 <- read_excel(
  "data/raw/mn_dnr_data/MN_IsotopeData_03182025.xlsx",
  sheet = "2024", 
  col_types = c(rep("text", 6), "date", "numeric", "numeric", rep("text", 6), rep("numeric", 7), "text")) %>% 
  mutate(across(where(is.character), ~na_if(., "NA"))) # replace "NA" with NA

# small fixes to mn_feathers_2025

# seven MN reference birds from washington county, MN need lat longs, use county centroid
mn_counties <- tigris::counties(state = "MN", cb = TRUE) %>% 
  st_transform(4152) # project to MN NAD83 HARN for st_centroid()
  
wash_cent <- mn_counties %>% 
  filter(NAME == "Washington") %>% 
  st_centroid() %>% 
  st_transform(4326) %>% # put back into lat/long
  mutate(
    x = st_coordinates(.)[, 1],  
    y = st_coordinates(.)[, 2]
  )

mn_feathers_2025 <- mn_feathers_2025 %>% 
  mutate(feather = case_when( # fix a few of 'feather' values
    feather == "Primary1" ~ "P1",
    feather == "P3/P4" ~ "P3_4",
    TRUE ~ feather)) %>% 
  mutate(latitude = case_when( # assign county centroid value to 7 MN reference samples from Washington County without lat/longs
    source == "reference" & is.na(latitude) ~ wash_cent$y,
    TRUE ~ latitude),
    longitude = case_when(
      source == "reference" & is.na(longitude) ~ wash_cent$x,
      TRUE ~ longitude)
  )

# combine mn data with rest of data and drop the two samples we didn't get d2h for
updated_feather_results_df <- rbind(feather_results_df, mn_feathers_2025) %>% 
  drop_na(d2h_vs_vsmow_original)

## save new versions of the results dataframes -- one with all feathers and one with just reference feathers
updated_feather_results_df %>% write_rds("data/processed/cleaned_feather_results/feather_results_df.rds")

## reference feathers dataframe

# need to add a column to indicate the newly received reference samples (WI 2022/23 and MN newly received data)
# these have not been sent to Jackson and incorporated into the known-origin datasets used to update calibration equations

# get list of ids from new reference samples
wi_new <- updated_feather_results_df %>% 
  filter(source == "reference" & batch_id == "batch4") %>% 
  select(cornell_sample_id)

mn_new <- mn_feathers_2025 %>% 
  filter(source == "reference") %>% 
  select(cornell_sample_id)

new_ref <- rbind(wi_new, mn_new) %>% pull()

updated_reference_df <- updated_feather_results_df %>% 
  filter(source == "reference") %>% 
  mutate(new_sample_2025 = case_when(
    cornell_sample_id %in% new_ref ~ 1,
    TRUE ~ 0
  ))

sum(updated_reference_df$new_sample_2025) # n = 186, yes

write_rds(updated_reference_df, "data/processed/cleaned_feather_results/reference_data_df.rds")

table(updated_reference_df$species)

# very basic data exploration here ----------------------------------------

# look for NA values across dataframe
updated_feather_results_df %>% 
  select(where(is.numeric)) %>%
  pivot_longer(
    cols = everything(),
    names_to = 'variable',
    values_to = 'value'
  ) %>%
  filter(is.na(value)) %>%
  print(n=Inf)

# look at histograms
updated_feather_results_df %>% 
  select(where(is.numeric)) %>%
  pivot_longer(
    cols = everything(),
    names_to = 'variable',
    values_to = 'value'
  ) %>%
  ggplot() +
  geom_histogram(aes(x = value)) +
  facet_wrap(~variable, scales = 'free')

updated_feather_results_df %>% 
  select(where(is.character)) %>%
  select(-cornell_sample_id, -town, -comments) %>%
  pivot_longer(
    cols = everything(),
    names_to = 'variable',
    values_to = 'value'
  ) %>%
  group_by(variable, value) %>%
  tally() %>%
  ggplot() +
  geom_col(aes(value, n)) +
  facet_wrap(~variable, scales = 'free') +
  coord_flip()
