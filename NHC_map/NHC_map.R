################################
# Make Map for hypoxia manuscript
# A Carter
# 5.20.2020


#install.packages("nhdplusTools") #uncomment and run!
library(nhdplusTools)
library(tidyverse)
library(sf)
library(tmap)
library(rgdal)
library(maps)
options(stringsAsFactors = FALSE)

setwd(hypox_projdir)
# Get sample station locations
sites <- read.csv("NHC_map/NC_synopticSamplingSites.csv", header=T, stringsAsFactors = F)
sites_sf <- sites[c(1:12),1:3] %>% 
  st_as_sf(coords=c("Long","Lat"),remove=F, crs=4326)
wwtp_sf <- sites[15,1:3] %>% 
  st_as_sf(coords=c("Long","Lat"),remove=F, crs=4326)

# site_MCconf <- NHC_sites[14,1:3] %>%
#   st_as_sf(coords=c("Long","Lat"), remove=F, crs=4326)


# NHC_gage_id <- "02097314" #USGS gage at blands
# cur_gage_id <- NHC_gage_id
# cur_id_string<- paste0("USGS-",cur_gage_id) #just add prefix to string
# cur_site <- list(featureSource = "nwissite",  #prepping below web quirie below, this just indicates the kind of data we have
#                  featureID = cur_id_string) 
# cur_line <- navigate_nldi(cur_site, "UT", "") #returns all upstream flowlines
# 
# #Flowlines above are just the ID's . Download everything!
# save_file_name <- "NHC_NHD_subset.gpkg" #gkpg is an open source geospatial format
# subset_nhdplus(cur_line$nhdplus_comid, 
#               save_file_name, 
#               "download")

cur_nhd <- st_read("NHC_map/NHC_NHD_subset.gpkg") 
 
 
# NHC <- sites[sites$site=="NHC",1:3]%>%
#   st_as_sf(coords=c("Long","Lat"),remove=F, crs=4326)

# attach comids to each sample station
comid_points<- rep(NA, 13)
for(i in 1:nrow(sites_sf)){
  comid_points[i]<- discover_nhdplus_id(sites_sf[i,])  
}


NHC_mainstem <- cur_nhd %>% 
  filter(gnis_name=="Mud Creek" | gnis_name=="New Hope Creek"|comid==8888400)

tm_shape(NHC_mainstem) + tm_lines()

UNHC_NHD <- NHC_mainstem %>%
  filter(hydroseq > 250019295 &gnis_name=="New Hope Creek")
UMC_NHD <- NHC_mainstem %>%
  filter(hydroseq > 250087220 & gnis_name=="Mud Creek" )

longitudinal_transect <- NHC_mainstem %>%
  filter(!comid %in% UNHC_NHD$comid)%>%
  filter(!comid %in% UMC_NHD$comid)

sample_site_NHD_reaches <- tibble(comid=comid_points) %>%
  left_join(NHC_m)%>% st_as_sf()

siteDat <- sample_site_NHD_reaches %>% 
  select(comid, streamorde, pathlength, slope, totdasqkm )
siteDat$site <- sites_sf$site

# read in NHC watershed
nhc_ws <- st_read("ws_shapefiles/nhc/wsboundary.shp",stringsAsFactors=FALSE)
mud_ws <- st_read("ws_shapefiles/mud/wsboundary.shp",stringsAsFactors=FALSE)

# # Get a basemap
# library(OpenStreetMap)
# om<-openmap(upperLeft = c(36.031,-79.1657), lowerRight = c(35.8819,-78.8658), type="esri-topo")
# omr <- raster::raster(om)
# plot(om)
tmap_mode("view")

tm_shape(nhc_ws)+tm_polygons(alpha=0, border.col="black", lwd=.5)+
  #  tm_shape(mud_ws)+tm_polygons(alpha=0, border.col="black",lwd=.5)+
  tm_shape(cur_nhd)+tm_lines(col = "grey60") +
  tm_shape(longitudinal_transect) + tm_lines(lwd=2)+
  tm_shape(sites_sf)+tm_dots(col="brown3", size=.05)+
  #  tm_shape(wwtp_sf)+tm_markers(shape=3, col="lightblue",size=.05)+
  tm_scale_bar(text.size = 1, position = "left") 

tmap_mode("plot")
par(bg=NA)
map<-tm_shape(nhc_ws)+tm_polygons(alpha=0, border.col="black", lwd=.5)+
#  tm_shape(mud_ws)+tm_polygons(alpha=0, border.col="black",lwd=.5)+
  tm_shape(cur_nhd)+tm_lines(col = "grey60") +
  tm_shape(longitudinal_transect) + tm_lines(lwd=2)+
  tm_shape(sites_sf)+tm_dots(col="brown3", size=.05)+
#  tm_shape(wwtp_sf)+tm_markers(shape=3, col="lightblue",size=.05)+
  tm_scale_bar(text.size = 1, position = "left") +
  tm_compass(type="arrow",position=c("right","bottom", show.labels=3))+
  tm_layout(frame=FALSE, bg.color="transparent")
tmap_save(map, filename="NHCmap_scalebar.eps", bg="transparent")


# Plot of longitudinal transect
long_sites_sf <- sites_sf[sites_sf$site!="MC751",]
tmap_mode("view")

tm_shape(cur_nhd)+tm_lines(col = "grey80") +
  tm_shape(longitudinal_transect) + tm_lines(lwd=2)+
  tm_shape(long_sites_sf)+tm_dots(col="brown3", size=.05)+
  #  tm_shape(wwtp_sf)+tm_markers(shape=3, col="lightblue",size=.05)+
  tm_scale_bar(text.size = 1, position = "left") 

# Plot North carolina with piedmont shape


pied <- readOGR(dsn="ncpiedmont_shape",
                layer="Piedmont_shape")

par(bg=NA)
png("NCmap.png",bg="transparent", type="windows")
map('state',region='North Carolina',fill=TRUE, col="white",bg="transparent",lwd=2)
  plot(pied,
     add=TRUE, 
     col="grey90")
  points(wwtp_sf$Long, wwtp_sf$Lat, col="brown3", pch=22, cex=3)
dev.off()