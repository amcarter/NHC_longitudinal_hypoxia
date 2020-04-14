# Pull data from streampulse website
# AM Carter
# 23 Sept 2019
# Adapted from:

#Basic StreamPULSE data processing pipeline
#Updated 10/29/18
#Contact Mike Vlah (vlahm13@gmail.com) with questions or comments.

# Update the streamPULSE package from github
# The StreamPULSE package is in development and changes frequently!
# If something doesn't work as expected, first try reinstalling.

#library(devtools)
#install_github('streampulse/StreamPULSE', dependencies=TRUE)

# Load packages.
library(StreamPULSE)
library(streamMetabolizer)
library(tidyr)
library(dygraphs)
library(xts)
library(lubridate)

# View all available sites and their metadata
# query_available_data(region='all')

# View all variables at a site, and full available time range for that site.
# Note that USGS depth and discharge data may be available for sites that
# have associated USGS gage IDs, even if depth and discharge do not appear among
# the variables returned here. If USGS data are available, they will be acquired
# automatically when you use prep_metabolism below. Likewise, air pressure and
# PAR estimates will be automatically acquired below, if necessary.
query_available_data(region='NC')

# Select site and date range for which to acquire StreamPULSE data.
# site_code is a combination of regionID and siteID
NHCsites2018 <- c('MC751','MC3','MC2','MC1','NHC5',
                  'NHC4','NHC3','NHC2','NHC1')
NHCsites2019 <- c('NHC','PM','CBP','WB','WBP','PWC','UNHC')

# Create Dataframe for site data - for now do this by generating the df for the first site
#  and renaming it, then use a loop to fill the rest.

# Download data from streampulse.
 start_date <- ymd_hms("2018-06-12 13:45:00 UTC")
 end_date <- ymd_hms("2018-07-06 02:45:00 UTC")

# Pull Air Pressure data from NHC core site
airP <- request_data(sitecode="NC_NHC", variables = "AirPres_kPa",
                    startdate=start_date, enddate=end_date)$data %>%
  spread(key = "variable", value="value") %>% select(c("DateTime_UTC", "AirPres_kPa"))
alldat <- airP

# Use this to convert water pressure data to water level data using sensor offset data:
# Note, as of 11/10/2019 this is a fake dataset that needs to be updated based on the field notes
sensorOffsets <- read.csv("data/sensorOffsets.csv", header=T)

# This for loop does a billion things:

for(site in NHCsites2018){
  site_code <- paste('NC',site,sep="_")
  
  # Get list of variables and start and end date for each site from the database
  mdat <- query_available_data(region = 'NC',site = site)
  variables <- as.vector(mdat$variables$variables)
   #end_date <- mdat$datebounds$lastRecord
  
  #Download data from SP
  dat <- request_data(sitecode=site_code, variables = variables,
                      startdate=start_date, enddate=end_date)
  w <- which(dat$data$flagtype== "Bad Data")
  dat$data$value[w]<-NA
  dat<- dat$data[,c(1,4,5)]%>% spread(key = "variable",value = "value")
  
  # attach air pressure data to calculate percent sat and depth
  tmp <- full_join(airP, dat, by = "DateTime_UTC")
  
  # Use function from Stream Metabolizer to calculate percent saturation
  dat$DO.sat <- calc_DO_sat(tmp$WaterTemp_C, tmp$AirPres_kPa)
  if("Level_m" %in% colnames(dat)){
    dat$Level_m <- dat$Level_m + sensorOffsets$offset.m[which(sensorOffsets$site==site)]
  }
  if("WaterPres_kPa" %in% colnames(dat)){
      dat$Level_m <- (tmp$WaterPres_kPa-tmp$AirPres_kPa)*0.10197 +# m of water per kPa
        sensorOffsets$offset.m[which(sensorOffsets$site==site)]
  }
  
  colnames(dat) <- c("DateTime_UTC",paste(site, colnames(dat[,-1]), sep = "."))
  
  alldat <- full_join(alldat,dat, by = "DateTime_UTC")
}


#Create Time series object of just the DO values to plot with DY graphs
DOdat <- data.frame(alldat[,c(1,2,
                              grep("*.DO_mgL",colnames(alldat)), 
                              grep("*.WaterTemp_C", colnames(alldat)),
                              grep("*.Level_m",colnames(alldat)),
                              grep("*.DO.sat", colnames(alldat)))])



write.csv(DOdat, file = "data/allsites_NHC2018DOdata.csv", row.names=F)

#DOdat.xts <- xts(alldat[,grep("*.DO_mgL", colnames(alldat))],order.by=alldat[,1])
#dygraph(DOdat.xts) %>% dyRangeSelector()



# Summarize by day
DOdat$date <- date(DOdat$DateTime_UTC)
DOsummarystats <- function(site, DOdat){
  daily <- DOdat[,c("date",paste0(site,".DO_mgL"))] %>% 
    rename("DO"=paste0(site,".DO_mgL"))%>%
    group_by(date) %>%                  
    summarize(mean =mean(DO,na.rm = T),
              min = min(DO, na.rm = T),
              max = max(DO,na.rm = T),
              amp = (max-min))
}
DOdaily <- DOdat %>%
  group_by(date)%>%
  summarize(
    for(i in (nrow(DOdat)-1)) {
      mean =mean(DO,na.rm = T),
      min = min(DO, na.rm = T),
      max = max(DO,na.rm = T),
      amp = (max-min)
    }      
  ) 

MC1daily <- DOsummarystats("MC1", DOdat)
MC2daily <- DOsummarystats("MC2", DOdat)
MC3daily <- DOsummarystats("MC3", DOdat)
MC751daily <- DOsummarystats("MC751", DOdat)
NHC1daily <- DOsummarystats("NHC1", DOdat)
NHC2daily <- DOsummarystats("NHC2", DOdat)
NHC3daily <- DOsummarystats("NHC3", DOdat)
NHC4daily <- DOsummarystats("NHC4", DOdat)
NHC5daily <- DOsummarystats("NHC5", DOdat)

startdate <- date("2018-06-12")
enddate <- date("2018-07-05")

startdate <- date("2018-10-02")
enddate <- date("2018-12-01")

plotDOstats <- function(MC1daily, site, startdate, enddate){
  par(mfrow = c(2,1))
  par(mar = c(0,5,0,0), oma = c(3,0,.5,.5))
  plot(MC1daily$date, MC1daily$mean, type = "l", cex = 1.2, ylab = "DO mgL",
       ylim = c(0,12), xaxt='n',xlim = c(startdate,enddate))
  mtext(site, side = 3, line = -1.2, adj = 0.03, cex = 1.2)
  lines(MC1daily$date, MC1daily$min, col = "grey60")
  lines(MC1daily$date, MC1daily$max, col = "grey60")
  plot(MC1daily$date, MC1daily$amp, type = "l", cex = 1.2,ylab = "DO amp",
       ylim = c(0,10), xlim = c(startdate,enddate), xaxt = 'n', yaxs="i")
  axis.Date(1, at=seq(startdate, enddate, by="weeks"), format="%m-%d")
  
  }
par(mfrow = c(4,2), cex = 1,
    mar = c(0,0,0,0), oma = c(3,3,.5,.5))

plotDO <- function(MC1daily, site, startdate,enddate)
  {
  plot(MC1daily$date, MC1daily$mean, type = "l", cex = 1.2, 
     axes=F, ylim = c(0,12), xlim = c(startdate,enddate))
  box(col="grey60")
  mtext(site, side = 3, line = -.8, adj = 0.03, cex = 0.6)
  lines(MC1daily$date, MC1daily$min, col = "grey60")
  lines(MC1daily$date, MC1daily$max, col = "grey60")
  axis(2, col = "grey40", col.axis = "grey20")
  plot(MC1daily$date, MC1daily$amp, type = "l", cex = 1.2, axes = F,
       ylim = c(0,12), xlim = c(startdate,enddate))
  box(col="grey60")
  }

#polygon(c(MC1daily$date, rev(MC1daily$date)), c(MC1daily$amp, rep(0, nrow(MC1daily))))
plotDO(MC751daily,"MC751", startdate, enddate)
plotDO(MC3daily,"MC3", startdate, enddate)
plotDO(MC2daily,"MC2", startdate, enddate)
plotDO(MC1daily, "MC1", startdate, enddate)
plotDO(NHC5daily,"NHC5", startdate, enddate)
plotDO(NHC4daily,"NHC4", startdate, enddate)
plotDO(NHC3daily,"NHC3", startdate, enddate)
plotDO(NHC2daily,"NHC2", startdate, enddate)
plot(NHC1daily$date, NHC1daily$mean, type = "l", cex = 1.2, 
     axes=F, ylim = c(0,12), xlim = c(startdate,enddate))
box(col="grey60")
mtext("NHC1", side = 3, line = -1, adj = 0.03, cex = 0.6)
lines(NHC1daily$date, NHC1daily$min, col = "grey60")
lines(NHC1daily$date, NHC1daily$max, col = "grey60")
axis(2, col = "grey40", col.axis = "grey20")
axis.Date(1, at=seq(startdate, enddate, by="weeks"), format="%m-%d")
plot(MC1daily$date, MC1daily$amp, type = "l", cex = 1.2, axes = F,
     ylim = c(0,12), xlim = c(startdate,enddate))
box(col="grey60")
axis.Date(1, at=seq(startdate, enddate, by="weeks"), format="%m-%d")


plotDO(NHC1daily,"NHC1")
axis(1, col = "grey40", col.axis = "grey20")
for (i in 1:6) {
  plot(1, axes = FALSE, type = "n")
  mtext(letters[i], side = 3, line = -1, adj = 0.9, cex = 0.6,
         col = "grey40")
  if (i %in% c(4, 5, 6))
    , at = seq(0.6, 1.2, 0.2))
  if (i %in% c(1, 4))
    , at = seq(0.6, 1.2, 0.2))
  box(col = "grey60")
 }
mtext("x axis", side = 1, outer = TRUE, cex = 0.7, line = 2.2,
         col = "grey20")
mtext("y axis", side = 2, outer = TRUE, cex = 0.7, line = 2.2,
         col = "grey20")

                  