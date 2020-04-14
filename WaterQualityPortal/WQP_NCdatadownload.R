###################################################################
# Water Quality Portal Data Download
# Dissolved Oxygen Data from North Carolina Piedmont

# A Carter
# Created: 2020-03-30

# install.packages("dataRetrieval")
library(dataRetrieval) # interact with NWIS and WQP data
library(tidyverse)
library(leaflet)
library(streamMetabolizer)
library(streamPULSE)
library(sf)
library(rgdal)
library(rgeos)


# Find outline of NC piedmont

dat <- readWQPdata(statecode="NC", characteristicName = c("Oxygen","Dissolved oxygen (DO)"), siteType="Stream", sampleMedia="Water")


DOdat <- dat %>% select(ActivityStartDate, ActivityStartTime.Time, ActivityStartTime.TimeZoneCode,
                  ActivityStartDateTime,MonitoringLocationIdentifier, CharacteristicName,
                  ResultMeasureValue, ResultMeasure.MeasureUnitCode)
####


# Remove data with unidentified units or the wrong sample media:


DOdat <- DOdat[which(DOdat$ResultMeasure.MeasureUnitCode=="mg/l"),]

w <- which(DOdat$ResultMeasureValue >20)
DOdat <- DOdat[-w,] %>% rename(DO_mgl=ResultMeasureValue) %>%
  select(-CharacteristicName, -ResultMeasure.MeasureUnitCode)

siteInfo <- attr(DOdat, "siteInfo") %>% select(MonitoringLocationIdentifier, latitude=dec_lat_va, longitude=dec_lon_va,
                                               CRS=HorizontalCoordinateReferenceSystemDatumName)

DOdat <- left_join(DOdat, siteInfo, by="MonitoringLocationIdentifier")
sites <- unique(DOdat$MonitoringLocationIdentifier)
startdate <- min(DOdat$ActivityStartDate, na.rm=T)

tempdat <- readWQPdata(siteid=sites, characteristicName="Temperature, water", startDateLo=startdate, siteType="Stream", sampleMedia="Water")
tempdat <- tempdat[which(tempdat$ResultMeasure.MeasureUnitCode=="deg C"),] %>% select(ActivityStartDate, ActivityStartTime.Time, MonitoringLocationIdentifier, temp_C=ResultMeasureValue)
tempdat$temp_C[which(tempdat$temp_C >39)]<- NA
alldat <- left_join(DOdat, tempdat, by=c("ActivityStartDate","ActivityStartTime.Time","MonitoringLocationIdentifier"))


#################################################
# Clip to NC piedmont using a shape file
pied_sf <- readOGR("Piedmont_shape.shp")
spatial.dat<- alldat
coordinates(spatial.dat)<- ~longitude + latitude
proj4string(spatial.dat) <- proj4string(pied_sf)

pied_points <- over(spatial.dat, pied_sf)
w<- which(!is.na(pied_points$count))

plot(pied_sf)
points(alldat$longitude, alldat$latitude)
points(alldat$longitude[w], alldat$latitude[w], col = "red")

dat_piedmont <- alldat[w,]
unique(dat_piedmont$MonitoringLocationIdentifier)

write_csv(dat_piedmont, "piedmontWQP_DOmgl.csv")
dat_piedmont <- read_csv("piedmontWQP_DOmgl.csv")


# convert all of the years
dat_piedmont$ActivityStartDateTime <- ymd_hms(paste0(as.character(dat_piedmont$ActivityStartDate), 
                                                     " ",as.character(dat_piedmont$ActivityStartTime.Time)), 
                                              tz="EST")
dat_piedmont$DateTime_UTC <-round_date(with_tz(dat_piedmont$ActivityStartDateTime,tz="UTC"),
                                       unit="15 minutes")

dat_piedmont <- dat_piedmont%>% select(-ActivityStartDate,-ActivityStartTime.Time,
                                       -ActivityStartTime.TimeZoneCode, -ActivityStartDateTime)
dat_piedmont <- dat_piedmont[-which(is.na(dat_piedmont$DateTime_UTC)),]
############################################################################
# get air pressure for all of the points in the piedmont
dat_piedmont$year <- year(dat_piedmont$DateTime_UTC)
dat_piedmont$airPres_kPa <- NA
dat_piedmont$air_temp_C <- NA

write_csv(dat_piedmont, "piedmontWQP_airP.csv")
saveRDS(dat_piedmont, "piedmontWQP_airP.rds")
dat_piedmont<- readRDS("piedmontWQP_airP.rds")
sites <- unique(dat_piedmont$MonitoringLocationIdentifier)
revisit = data.frame()

for(i in 1:length(sites)){
  w <- which(dat_piedmont$MonitoringLocationIdentifier==sites[i])
  tmp <- dat_piedmont[w,]
  yrs <- unique(tmp$year)
  airP <- data.frame(DateTime_UTC=ymd_hms("2020/01/01 00:00:00"), air_kPa = 1, air_temp=1)
  for(j in 1:length(yrs)){
    advance_loop<- FALSE
    tryCatch({
      airP_y<-StreamPULSE:::FindandCollect_airpres(tmp$latitude[1], tmp$longitude[1],
                                               ymd_hms(paste0(yrs[j],"/01/01 00:15:00")), 
                                               ymd_hms(paste0(yrs[j],"/12/31 23:45:00")))
      
    },error=function(e){
      revisit <<-bind_rows(revisit, data.frame(i,j))
      advance_loop <<-TRUE
    })
    if(advance_loop) next 
    airP <- rbind(airP, airP_y)
    }

  airP$MonitoringLocationIdentifier<- sites[i]
  tmp <- left_join(tmp, airP, by=c("MonitoringLocationIdentifier", "DateTime_UTC"))
  dat_piedmont$air_temp_C[w]<- tmp$air_temp 
  dat_piedmont$airPres_kPa[w]<- tmp$air_kPa
  print(paste0(i/length(sites)*100, "% done"))
  if(i %% 50 == 0){saveRDS(dat_piedmont, "piedmontWQP_airP.rds")}
  }
t <- sites[1:85]
tmp <- dat_piedmont[dat_piedmont$MonitoringLocationIdentifier %in% t,]
length(which(is.na(tmp$airPres_kPa)))/nrow(tmp)
plot(tmp$temp_C, tmp$air_temp_C)


###########################################################
# Convert DO data to percent saturation using StreamMetabolizer package
# Pressure must be converted to millibars first: 1 kPa = 10 mbar
dat_piedmont$DO_sat<- calc_DO_sat(dat_piedmont$temp_C, 10*dat_piedmont$airPres_kPa)
dat_piedmont$DO_persat <- dat_piedmont$DO_mgl/dat_piedmont$DO_sat 
dat_piedmont$month <- months(dat_piedmont$DateTime_UTC)

plot(dat_piedmont$DO_mgl, dat_piedmont$DO_persat)

annual <- dat_piedmont %>%group_by(year)%>%
  summarize(DO_persat.mean=mean(DO_persat, na.rm=T),
            DO_persat.sd=sd(DO_persat, na.rm=T))
plot(annual$year, annual$DO_persat.mean, ylim = c(0,1.2), type="l")
polygon(x=c(annual$year, rev(annual$year)), 
          y=c(annual$DO_persat.mean+annual$DO_persat.sd, rev(annual$DO_persat.mean-annual$DO_persat.sd)),
            col=alpha("black", .2))

hist(dat_piedmont$DO_persat, xlim = c(0,1.5), ylim = c(0,3100))
par(new=T)
hist(dat_piedmont$DO_persat[dat_piedmont$month=="October"], xlim = c(0,1.5), 
     ylim = c(0,3100), col = "brown3", axes=F)

plot(dat_piedmont$DateTime_UTC, dat_piedmont$DO_persat, col = alpha("black",0.2), pch=19, ylim = c(0,1.5))
points(dat_piedmont$DateTime_UTC[dat_piedmont$month=="October"],
       dat_piedmont$DO_persat[dat_piedmont$month=="October"],pch=19, col=alpha("brown3",.6) )
plot(annual$year, annual$DO_persat.mean, type="l", ylim = c(0,1.5))

ncSummary <- DOdat %>%
  group_by(MonitoringLocationIdentifier) %>%
  summarise(count=n(),
            max = max(value, na.rm = TRUE)) %>%
  filter(count > 10) %>%
  arrange(-count) %>%
  left_join(siteInfo, by = "MonitoringLocationIdentifier")


col_types <- c("darkblue","dodgerblue","green4","gold1","orange","brown","red")
leg_vals <- unique(as.numeric(quantile(ncSummary$max,
                                       probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
pal = colorBin(col_types, ncSummary$max, bins = leg_vals)
rad <-3*seq(1,4,length.out = 16)
ncSummary$sizes <- rad[as.numeric(cut(ncSummary$count, breaks=16))]

leaflet(data=ncSummary) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(~dec_lon_va,~dec_lat_va,
                   fillColor = ~pal(max),
                   radius = ~sizes,
                   fillOpacity = 0.8, opacity = 0.8,stroke=FALSE,
                   popup=~station_nm) %>%
  addLegend(position = 'bottomleft',
            pal=pal,
            values=~max,
            opacity = 0.8,
            labFormat = labelFormat(digits = 1),
            title = 'Max Value')

