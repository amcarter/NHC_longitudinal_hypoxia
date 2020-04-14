# Delineate MC and NHC watersheds, determine watershed characteristics

# watershed delineation with streamstats package ####
# Code from mvlah

# library(devtools)
# install_github('markwh/streamstats')
library(streamstats)
library(rgdal)

sites <- read.csv("NC_synopticSamplingSites.csv", header=T)
MCout <- sites[sites$site=="MCconf",2:3]
NHCout <- sites[sites$site=="Blands",2:3]

MC = streamstats::delineateWatershed(MCout$Long, MCout$Lat,crs=4326,
                                     includeparameters="true",includeflowtypes = "true")
streamstats::writeShapefile(MC,layer = "MC_bound",  what="boundary")

#save boundary as geojson, read back in as spatial dataframe, write as shapefile.
#let me know if you find a less roundabout way to do this.
 tpf = tempfile(fileext='.geojson')
streamstats::writeGeoJSON(MC, file='MC.geojson', what='boundary')
spatialdf = rgdal::readOGR("MC")
unlink(tpf) #remove tempfile (unnecessary)
rgdal::writeOGR(spatialdf, dsn='path/where/shapefile/components/will/go',
    layer='shapefile_name_no_extension', driver='ESRI Shapefile')

#create interactive map of ws boundary
#seems like saveWidget only saves to cwd. precede with setwd() to save to subdir
webfile = streamstats::leafletWatershed(ws_bound)
htmlwidgets::saveWidget(webfile, selfcontained=FALSE,
    file='file.html')
