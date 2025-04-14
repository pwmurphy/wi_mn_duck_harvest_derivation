
library(assignR)
library(isocat)
library(raster)
library(terra)
library(sf)
library(sp)
library(rasterVis)
library(rnaturalearth)
library(dplyr)
library(ggplot2)
library(mapview)
library(rasterVis)
library(viridis)

library(readr)
set.seed(2189)

## WDNR Isocat Assignment Demo
## Aug 27, 2024

# Isotope dataframe with individual VSMOW values (and metadata)
#d2Hf <- readRDS("C:/Users/SStemaly/OneDrive - LSU AgCenter/LA_Waterfowl_Origins/Other_Outputs/2023_FeatherSamples/2023_IsotopeSamples.rds") 
d2Hf <- read_rds("output/reference_data_df.rds") 
# Filtering to one species
GADW <- d2Hf %>%
  filter(Species == 'GADW') 

WODU <- d2Hf %>%
  filter(species == 'WODU') 

# Isoscapes. We used mean annual 
isoscape <- rast("data/assignment_data/precipitation_isoscape_MA/d2h_MA.tif") 
isoscape.se <- rast("data/assignment_data/precipitation_isoscape_MA/d2h_se_MA.tif")

# Calibration of isoscape. Make sure calibration equation matches the isoscape selected.
cal.isoscape <- isoscape * 0.6 - 71.7 # MA dabblers equation - From Kusack et al. (2023)
sd.resid <- 17.2  

# Species breeding range
WODU.Breeding.Range.rds <- read_rds(file = "data/assignment_data/spp_assignment_ranges/mall_bli_assignment_range.rds") %>% 
  st_as_sf() %>% 
  st_transform(crs(cal.isoscape))

# Set breeding range crs to the isoscape
# desired_crs <- crs(cal.isoscape)
# GADW.Breeding.Range.rds <- st_transform(GADW.Breeding.Range.rds, crs = desired_crs)
WODU.Breeding.Range.rds <- vect(WODU.Breeding.Range.rds)

# Crop the isoscape to the breeding range
cal.isoscape <- mask(cal.isoscape, WODU.Breeding.Range.rds) %>% 
  crop(ext(WODU.Breeding.Range.rds)) 

cal.isoscape.se <- mask(isoscape.se, WODU.Breeding.Range.rds) %>% 
  crop(ext(WODU.Breeding.Range.rds)) 

# Assignment
prob <- isotopeAssignmentModel(ID = WODU$cornell_sample_id, 
                               isotopeValue = WODU$d2h_vs_vsmow_original, 
                               SD_indv = sd.resid, 
                               precip_raster = raster(cal.isoscape), 
                               precip_SD_raster = raster(cal.isoscape.se)) #Output is a large rasterstack
plot(prob[[1]])


###### PM code #########
#####Binary assignment with 2:1 odds ratio ####
odds <- 0.33333333

d = na.omit(data.frame(values(s)))
e <- reshape(d, direction = "long", v.names = "density", timevar = "layer", varying = c(1:ncol(d)), sep = "") # reshape the data frame, creating a variable called "layer" to act as a factor

head(e)

fun <- function(x) predict(smooth.spline(cumsum(sort(x)), sort(x), spar = 0.1), odds) # Function to fit splines to estimate the probability densities at which the threshold cumulative probability is reached

cutoffs <- by(e$density, e$layer, fun) # applies function to fit splines for each layer (bird)
cutoffs <- as.matrix(unlist(cutoffs))
cutoffs <- as.matrix(cutoffs[-which(cutoffs[,1] == odds)])
colnames(cutoffs) <- "cutoff"

# view one bird
temp <- data.frame(cumsum(sort(e$density[e$layer == 3])),sort(e$density[e$layer == 3]))
names(temp) <- c("cumsum","prob")

plot(smooth.spline(temp$cumsum, temp$prob, spar = 0.1), ylab = "Posterior Probability (Ascending)", xlab = "Cumulative Posterior Probability", cex = 0.1)
abline(v = odds, lty = 2, col = "red")
abline(h = cutoffs[3], lty = 2, col = "blue")
text(0.15, max(temp$prob) - (max(temp$prob) - min(temp$prob))/20, paste("Cutoff: ", round(cutoffs[3], 8), sep = ""))





# Run the clusters based on the prob surfaces
myMatrix <- simmatrixMaker(prob)
cS <- clusterSimmatrix(simmatrix = myMatrix, dist_mthd = "correlation", hclust_mthd = "average")

# Determine the clusters based on the au probabilities. Cutting at 0.5 based on Kusack et al. 2022
cut.prob <- pvclust::pvpick(cS, alpha = 0.5, pv = "au")

# Create a new dataframe from our original dataframe (GADW) and assign the clusters based on the cut.prob object
# I used this approach to assign the clusters, rather than merging dataframes, because the cut.prob object is a list
# There are probably better ways to combine these data
# Finally, merge the clusters 
GADW$cluster <- NA
for(g in seq_len(length(cut.prob$clusters))) { GADW$cluster[GADW$cornell_sample_id %in% cut.prob$clusters[[g]]] <- g } #This assigns each individual to a cluster

# Create binary regions
prob <- rast(prob) # convert prob to a SpatRaster

# Despite Isocat, qtlRaster is still the best function to get binary surfaces
s <- assignR::qtlRaster(prob, threshold = 0.66, thresholdType = "prob") # this creates maps of the top 66% likelihood assignment for each individual
s <- s * 1  

bin.surface <- s %>% project("EPSG:5070", method = 'near') # Always show my final depictions in a conic crs -- SS used the Albers Equal Area Conic projection
# When reprojecting the binary surfaces, always use nearest neighbor interpolation, otherwise you get values in between 0-1 at the edges

# Group by cluster
# Here the origins object will have a separate layer for each cluster, where all the individuals are summed within in cluster
origins <- list()
for(g in 1:length(unique(GADW$cluster[is.na(GADW$cluster) == F]))) { # Note: sometimes samples are not assigned to any cluster, hence the is.na()
  clustStack <- subset(bin.surface, GADW$cornell_sample_id[GADW$cluster == g & is.na(GADW$cluster) == F])
  origins[[g]] <- app(clustStack, fun=function(k){sum(k)})
}

origins <- rast(origins)
names(origins) <- paste("Cluster", 1:max(GADW$cluster[is.na(GADW$cluster) == F]), sep = ".")

# If any samples are too different to cluster, we can depict them separately
if(sum(is.na(GADW$cluster)) >= 1) {
  origins<- c(origins, subset(bin.surface, GADW$cornell_sample_id[is.na(GADW$cluster) == T]))
}

# Finally I also make a population version, where all birds are combined
population.origins <- app(origins, fun=function(k){sum(k)})

# Raster with all the clusters, non-clustered birds, and the total surface 
origins.total <- c(origins, population.origins)
names(origins.total)[nlyr(origins.total)] <- "Total"

# ----------------------------------------------------------------------------------
# Prep for plotting

# Read in north america - making all CRS systems the same
northamerica.states <- ne_states(country =  c("Canada","United States of America"), returnclass = "sf") %>% 
  filter(!(name == "Hawaii")) %>%
  st_transform("+init=epsg:4326 +proj=longlat")

# Crop the merged boundary to the bounding box
bbox_coords <- st_bbox(c(xmin = -175, ymin = 5, xmax = -50, ymax = 83), crs = st_crs(northamerica.states))
cropped_NA <- st_intersection(northamerica.states, st_as_sfc(bbox_coords))
cropped_NA_trans <- st_transform(cropped_NA, crs = "+init=epsg:5070") # Clunky, but changing the crs to conic here
mapview(cropped_NA_trans)

Final.Northamerica <- as(cropped_NA_trans, "Spatial") # Changing class to sp object large spdf for plot layering later
mapview(Final.Northamerica)

# Convert GADW breeding range to match the population origin and north america CRS system
GADW.Breeding.Range.trans<- readRDS(file = "Other_Outputs/Species_Range_RDs_Files/GADW.Breeding.Range.rds") %>% 
  st_as_sf() %>%
  st_make_valid() %>%
  st_transform(crs = 5070) %>%
  summarise()

GADW.Breeding.Range.trans <- as(GADW.Breeding.Range.trans, "Spatial") # Changing class to sp object large spdf for plot layering later
mapview(GADW.Breeding.Range.trans)

## ------------ Combined origins for all included individuals
p.orig <- levelplot(population.origins, col.regions = c("#e6e6e6",viridis(n = 16, direction = -1)), 
                    margin = FALSE, scales = list(draw = FALSE), xlab = NULL, ylab = NULL, 
                    colorkey = list(title = "Count", title.gpar = list(col = "white"), title.control = list(side = "top"))) +
  latticeExtra::layer(sp.polygons(Final.Northamerica, size = 0.05)) +    # Layering on the spdfs we created
  latticeExtra::layer(sp.polygons(GADW.Breeding.Range.trans, col = "red", size = 0.01, lwd = 0.6))

p.orig

## ------------ Clustered origins

# Include sample sizes in plot
# Create a list that stores the sample size by cluster
samples.list <- list()
for (j in 1:length(unique(GADW$cluster[is.na(GADW$cluster) == F]))) { 
  samples.list[[j]] <- paste("n = ", table(GADW$cluster)[[j]], sep = "")
}

# Add n = 1 for any non-clustering birds
if(sum(is.na(GADW$cluster)) >= 1) {
  for (j in (length(samples.list) + 1):(length(samples.list) + sum(is.na(GADW$cluster)))) {
    samples.list[[j]] <- "n = 1"
  }
}

# Combined origins for each cluster
p.clust <- levelplot(origins, col.regions = c("#e6e6e6",viridis(n = 16, direction = -1)), 
                     margin = FALSE, scales = list(draw = FALSE), xlab = NULL, ylab = NULL, 
                     names.attr = paste(gsub(".", " ", names(origins), fixed=TRUE), " (", unlist(samples.list), ")", sep = ""),
                     par.settings = list(layout.heights = list(top.padding = 0, bottom.padding = 0),
                                         layout.widths = list(right.padding = 0, left.padding = 0),
                                         strip.background = list(col = "white")),
                     colorkey = list(title = "Count", title.gpar = list(col = "white"), title.control = list(side = "top"))) + # For some reason, the code is giving 2 titles, so I am making of them white
  latticeExtra::layer(sp.polygons(Final.Northamerica, size = 0.05)) +    # Layering on the spdfs we created
  latticeExtra::layer(sp.polygons(GADW.Breeding.Range.trans, col = "red", size = 0.01, lwd = 0.6))

p.clust


# ---------------------------------------------------------------------------------
# Mean Surface Cluster Map

# convert "prob" from a SpatRaster to RasterStack
probabilty_surface <- stack(prob)

# create aggregate surface of MEAN within-group probability of origin
# The below creates a cluster surface based on the MEAN, therefore there is some variance 
# associated with the boundaries of the cluster and we see that in the summary surfaces of all individuals within each cluster
meanSurfaces_all <- meanAggregateClusterProbability(
  indivIDs = GADW$cornell_sample_id, 
  clusters = GADW$cluster, 
  surfaces = probabilty_surface, 
  nClust = FALSE 
)

# Sets up thematic format for plot
gglayers <-  list(
  geom_tile(aes(fill = value)),
  coord_equal(),
  theme_bw(),
  scale_x_continuous(name = "Long", expand = c(0,0)),
  scale_y_continuous(name = "Lat", expand = c(0,0)),
  labs(title= "Cluster Classification")
)

# Specifies additional formatting 
ggProb <- list(
  facet_wrap(~ variable),
  scale_fill_gradient(name = "Probability\nOf Origin", low = 'darkblue', high = 'yellow')
)

# Actual plotting saved as an object with the addition of the formatting
cluster_layers <- gplot(meanSurfaces_all) + gglayers + ggProb

# Create a summary surface showing which RasterLayer in a Stack has the highest value at a given location
summary_map_all <- projectSummaryMaxSurface(surfaces = meanSurfaces_all, nClust = FALSE)
raster::plot(summary_map_all)

summary_map_all_projected <- projectRaster(summary_map_all, crs=crs(Final.Northamerica))

summary_map_all_projected_rast <- as(summary_map_all_projected, "SpatRaster")

# Plotting of the cluster surfaces
colors <- c("#A5DB36", "#2A788E")

plot(summary_map_all_projected, col = colors, main = "GADW Cluster Geographies",
     xlab = "longitude",
     ylab = "latitude")
plot(Final.Northamerica, add = T)

