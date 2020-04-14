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
NHCsites2018 <- c('Mud','MC751','MC3','MC2','MC1','NHC5',
                  'NHC4','NHC3','NHC2','NHC1','NHC','UNHC')
NHCsites2019 <- c('NHC','PM','CBP','WB','WBP','PWC','UNHC')

# Create Dataframe for site data - for now do this by generating the df for the first site
#  and renaming it, then use a loop to fill the rest.

# Download data from streampulse.
# These are the dates for which we have data at all sites
 start_date <- ymd_hms("2018-06-12 13:45:00 UTC")
 end_date <- ymd_hms("2018-07-06 02:45:00 UTC")

# Pull Air Pressure data from NHC core site
airP <- request_data(sitecode="NC_NHC", variables = "AirPres_kPa",
                    startdate=start_date, enddate=end_date)$data %>%
  spread(key = "variable", value="value") %>% select(c("DateTime_UTC", "AirPres_kPa"))

# Use this to convert water pressure data to water level data using sensor offset data:
# Note, as of 11/10/2019 this is a fake dataset that needs to be updated based on the field notes
sensorOffsets <- read.csv("data/sensorOffsets.csv", header=T)

# List of all the variables I might want from a given site:
variables <- c("DO_mgL", "satDO_mgL","Level_m", "SpecCond_uScm", "WaterPres_kPa", "WaterTemp_C")

getSPdata <- function(site, start_date, end_date, variables, airP, sensorOffsets){
  site_code <- paste('NC',site,sep='_')
  dat <- request_data(site_code, start_date, end_date, variables)
  # remove bad data
    w <- which(dat$data$flagtype== "Bad Data")
    dat$data$value[w]<-NA
  dat<- dat$data[,c(1,4,5)]%>% spread(key = "variable",value = "value")
  # attach air pressure data to calculate percent sat and depth
  dat <- full_join(airP, dat, by = "DateTime_UTC")

  # Use function from Stream Metabolizer to calculate percent saturation
  if(!("satDO_mgL" %in% colnames(dat))){
    dat$satDO_mgL <- calc_DO_sat(dat$WaterTemp_C, dat$AirPres_kPa*10)
    dat$persatDO <- dat$DO_mgL/dat$satDO_mgL
  }
  
  # Determine if level data is already available, if not calculate it from
  # Water Pressure. Then add sensor offset values
  if("WaterPres_kPa" %in% colnames(dat)){
    dat$Level_m <- (dat$WaterPres_kPa-dat$AirPres_kPa)*0.10197 +# m of water per kPa
      sensorOffsets$offset.m[which(sensorOffsets$site==site)]
  }

  
  return(dat)
}
  
# create a dataframe for all fo the data together:
alldat <- data.frame()

for(site in NHCsites2018){
  dat <- getSPdata(site, start_date, end_date, variables, airP, sensorOffsets)
  write.csv(dat, file = paste0("data/raw/",site,".csv"),row.names = F)
  dat$site <- rep(site, nrow(dat))
  alldat <- bind_rows(alldat, dat)
}

write.csv(alldat, file = "data/raw/allSites_20180706.csv", row.names = F)




#Create Time series object of just the DO values to plot with DY graphs
# DOdat <- data.frame(alldat[,c(1,2,
#                               grep("*.DO_mgL",colnames(alldat)), 
#                               grep("*.WaterTemp_C", colnames(alldat)),
#                               grep("*.Level_m",colnames(alldat)),
#                               grep("*.DO.sat", colnames(alldat)))])
# 
# 
# 
# write.csv(DOdat, file = "data/allsites_NHC2018DOdata.csv", row.names=F)
# 
# #DOdat.xts <- xts(alldat[,grep("*.DO_mgL", colnames(alldat))],order.by=alldat[,1])
#dygraph(DOdat.xts) %>% dyRangeSelector()

