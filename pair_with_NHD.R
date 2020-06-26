###########################################################
#
# link longitudinal sampling sites to NHD flowlines


#install.packages("nhdplusTools") #uncomment and run!
library(nhdplusTools)
library(sf)
library(tidyverse)

setwd(hypox_projdir)

dat <- read_csv("2018 NHC MC Longitudinal/Longitudinal Survey May 2018/Longitudinal Data/LongitudinalSamples_2018May.csv")
dat <- dat[!(dat$streamSection %in% c("MCU_ConfAmVi", "MC_751")),]


long_sites <- dat %>% select(streamSection, Date, Time, distance_m, Latitude, Longitude, width_m, depth_m, velocity_ms, 
                             temp_C, DO_pctsat, DO_mgL, Cl.mgL, SO4.mgL, NO3.N.mgL, NH4.N.mgL, Habitat,
                             RoadCrossing ,WWTP, SampleStation)
long_sites$Longitude <- -long_sites$Longitude

long_sites_sf <- long_sites%>%
    st_as_sf(coords=c("Longitude","Latitude"),remove=F, crs=4326)

long_sites$comid <- NA

cur_nhd <- st_read("NHC_map/NHC_NHD_subset.gpkg") 
tm_shape(cur_nhd)+tm_lines(col = "grey60") +
  tm_shape(long_sites_sf)+tm_dots(col="brown3", size=.05)


for(i in 1:nrow(long_sites)){
  long_sites$comid[i]<- discover_nhdplus_id(long_sites_sf[i,])  
}


long_sites <- left_join(long_sites, cur_nhd[,c("comid","slope","streamorde")])
plot(long_sites$distance_m, long_sites$slope)
par(new=T)
plot(long_sites$distance_m, long_sites$DO_pctsat, type="l")


long_sites$Habitat[1:12] <- "Ri"
boxplot(slope~Habitat, data=long_sites)

write_csv(long_sites, "data/long_survey_withNHD.csv")

#########################################################################
# read in WQP data
WQPdat <- readRDS("WaterQualityPortal/piedmontWQP_airP_done.rds")
WQPdat_sf <- WQPdat %>% 
  st_as_sf(coords=c("longitude","latitude"),remove=F, crs=4326)

WQPdat$comid<- NA
for(i in 1:nrow(WQPdat_sf)){
  
  if(i %% 100 == 0) print(i/nrow(WQPdat))
  
  out = try(discover_nhdplus_id(WQPdat_sf[i,]))
  if(! 'try-error' %in% class(out)){
    WQPdat$comid[i]<- out
  }
}

write_csv(WQPdat)

comids <- unique(WQPdat$comid)
comids <- comids[-1183]
save_file_name <- "WaterQualityPortal/WQP_NHD_subset.gpkg" #gkpg is an open source geospatial format
subset_nhdplus(comids, 
               save_file_name, 
               "download")
all_nhd <- st_read(save_file_name)

dat <- all_nhd %>%
  select(comid, slope,hydroseq, streamorde)
WQPdat_nhd <- left_join(WQPdat, dat, by="comid")
WQPdat_nhd$slope[WQPdat_nhd$slope< -10]<- NA
boxplot(slope~streamorde, data=WQPdat_nhd)


WQPdat_flat <- WQPdat_nhd[WQPdat_nhd$slope<0.01,]
plot(WQPdat_flat$slope, WQPdat_flat$DO_mgl)
points(long_sites$slope, long_sites$DO_mgL, col=3, pch=19)
mm<- glm(DO_mgl ~slope, data=WQPdat_flat)
abline(7.89, 128, col=2)
mm <- glm(DO_mgl ~ slope + areasqkm, data=WQPdat_nhd)
summary(mm)


write_csv(WQPdat_nhd, "WaterQualityPortal/WQP_NHDlinked.csv")
