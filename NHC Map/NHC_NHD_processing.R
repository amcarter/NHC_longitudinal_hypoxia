######################
# Created by: Nick Bruns
#             10/29/2019
#
# Script for extracting NHD flowlines in NHC for AMC
# Consult these documents to understand NHD column values, their units, etc.
#  https://s3.amazonaws.com/nhdplus/NHDPlusV21/Data/NationalData/0Release_Notes_NationalData_Seamless_GeoDatabase.pdf


# Ok! First, we'll do the processing itself, triming NHD using the Jordan Lake HUC10 from the National Watershed Boundary Database.

# recquires these packages
 library(nhdplusTools)
 library(tidyverse)
 library(sf)
 library(tmap)

# Load GPS points of sensor stations
sites <- read.csv("NC_synopticSamplingSites.csv", header=T)
blands <- sites[sites$site=="Blands",1:3] %>%
  st_as_sf(coords=c("Long","Lat"),crs=4326)

sites<- sites[sites$site!=c("Blands","MCconf"),1:3]
glimpse(sites)
sites_sf <- sites %>% 
  st_as_sf(coords=c("Long","Lat"), crs=4326)


#another way to enforece brewing up data only when lacking!
# if(!file.exists("../data/in/jordan_lake_NHD.gpkg")){ 
LOAD_SAVED_DATA <- TRUE
if(!LOAD_SAVED_DATA){
  #load and extract HUC10 polygon
  st_layers(wbd_shape_dir)
  huc10_shapes <-st_read(wbd_shape_dir,"WBDHU10")
  
  jordan_lake_huc10_name <- "B Everett Jordan Lake-New Hope River" # found by filtering by state, and skimming the HUC10 NAME field
  jordan_lake_huc10_shape <- huc10_shapes %>% filter(NAME==jordan_lake_huc10_name) 
  
  #load the NHD, and use the polygon to extract it
  nhdplus_path("~/Documents/data_sets_larger/NHDPlusNationalData/NHDPlusV21_National_Seamless.gpkg") # takes a while the first time, as it slices up your geo-package 
  #and saves processed versions to your machine.
  # afte 1 time, it just returns a directory path.
  
  nhd_paths <- stage_national_data()
  network <- readRDS(nhd_paths$flowline)
  
  jordan_lake_NHD  <- network[jordan_lake_huc10_shape,] #this spatial operation takes a few minutes (5 or so)
  
  st_write(jordan_lake_NHD,"jordan_lake_NHD.gpkg")
  st_write(NHC_mainstem_NHD,"new_hope_creek_mainsteam_NHD.gpkg")
  st_write(jordan_lake_huc10_shape,"jordan_lake_huc10_shape.gpkg")
  
} else{
  jordan_lake_huc10_shape <- st_read("jordan_lake_huc10_shape.gpkg")
  NHC_mainstem_NHD <- st_read("new_hope_creek_mainsteam_NHD.gpkg")
  jordan_lake_NHD <- st_read("jordan_lake_NHD.gpkg")
}

glimpse(jordan_lake_NHD)

# Take a look at the full region NHD

tmap_mode("view")

jordan_lake_NHD %>% 
  # filter(StreamOrde>=2) %>%  
  mutate(StreamOrde=as.factor(StreamOrde)) %>% 
  tm_shape() + tm_lines(col="StreamOrde")  + 
  tm_shape(jordan_lake_huc10_shape) + tm_polygons(alpha=0)

# Extract the New Hope Creek mainstem flowlins

NHC_mainstem_NHD <- jordan_lake_NHD %>% 
  filter(GNIS_NAME=="New Hope Creek") 

# Use the field Path length (km's to outlet) to get position along mainstem

tm_shape(NHC_mainstem_NHD) + tm_lines(col='Hydroseq',style="cont")


# Now do the same for Mud Creek:
NHC_MC_mainstem_NHD <- jordan_lake_NHD %>%
  filter(GNIS_NAME=="Mud Creek" | GNIS_NAME=="New Hope Creek"|COMID==8888400)

MC_mainstem_NHD <- jordan_lake_NHD %>%
  filter(GNIS_NAME=="Mud Creek")

NHC_conf_COMID <- 8893236
NHC_conf_hydroseq <- jordan_lake_NHD$Hydroseq[which(jordan_lake_NHD$COMID==NHC_conf_COMID)]

UNHC_NHD <- jordan_lake_NHD %>%
  filter(Hydroseq > NHC_conf_hydroseq & GNIS_NAME=="New Hope Creek" )

MC_conf_COMID <- 8888406
MC_conf_hydroseq <- jordan_lake_NHD$Hydroseq[which(jordan_lake_NHD$COMID==MC_conf_COMID)]

UMC_NHD <- jordan_lake_NHD %>%
  filter(Hydroseq > MC_conf_hydroseq & GNIS_NAME=="Mud Creek" )

Mtrib <- jordan_lake_NHD %>%
  filter(COMID==8888400)

# Link sites to the nearest NHD reach
sites_sf %>%
  st_nearest_feature(NHC_MC_mainstem_NHD)

blands %>%
  st_nearest_feature(jordan_lake_NHD)

longitudinal_transect <- NHC_MC_mainstem_NHD %>%
  filter(!COMID %in% UNHC_NHD$COMID)%>%
  filter(!COMID %in% UMC_NHD$COMID)
tmap_mode("view")

tm_shape(longitudinal_transect) + 
  tm_lines(lwd=2) +
tm_shape(UNHC_NHD) +
  tm_lines(col="grey70",lwd=2)+
tm_shape(UMC_NHD)+
  tm_lines(col="grey70", lwd=2)+
tm_shape(sites_sf) +
  tm_dots(col="red", size=.04)+
tm_scale_bar(size = 1, position = "left") 

  #tm_dots(col="red",popup.vars=c("Site"))

export::graph2ppt("NHCtestMap")
# Inspect elevation and slope NHD values! 
# Probably preffered if plotted on same axess, here they are not.
# Again, use these docs to understand units, column values etc

# Note: if you give me lat-longs of your sites, we can link them to NHD reaches 
#       and use those to put vertical bars on the below elevation plots. 
#       We could do the same for the sample points.

NHC_MC_mainstem_NHD %>% 
  ggplot() + geom_line(aes(x=Pathlength,y=MINELEVSMO),col="red") +
  geom_line(aes(x=Pathlength,y=MAXELEVSMO),col="blue") +
  xlab("km to outlet") +
  xlab("distance to outlet (km)") +
  ggtitle("Downstream reach elevation. Blue=max, red=min")

NHC_MC_mainstem_NHD %>% 
  ggplot() + geom_line(aes(x=Pathlength,y=SLOPE)) +
  xlab("km to outlet")
