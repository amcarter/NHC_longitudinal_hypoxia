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
source("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/Code/Tools/timeseries_tools.R")

# View all available sites and their metadata
# query_available_data(region='all')

# View all variables at a site, and full available time range for that site.
# Note that USGS depth and discharge data may be available for sites that
# have associated USGS gage IDs, even if depth and discharge do not appear among
# the variables returned here. If USGS data are available, they will be acquired
# automatically when you use prep_metabolism below. Likewise, air pressure and
# PAR estimates will be automatically acquired below, if necessary.
# query_available_data(region='NC')

# Select site and date range for which to acquire StreamPULSE data.
# site_code is a combination of regionID and siteID

SPsites <- c("Mud","MC751","UNHC","NHC")
# Create Dataframe for site data - for now do this by generating the df for the first site
#  and renaming it, then use a loop to fill the rest.

# Download data from streampulse.
####################################################################
# Intersite data
# This paper will do intersite comparison of 2019 data and interannual at NHC from 2017-2019
 start_date <- ymd_hms("2019-01-01 05:00:00 UTC") # for intersite
 end_date <- ymd_hms("2020-01-01 05:00:00 UTC")

# Pull Air Pressure data from NHC core site
airP <- request_data(sitecode="NC_NHC", variables = "AirPres_kPa",
                    startdate=start_date, enddate=end_date)$data %>%
  spread(key = "variable", value="value") %>% select(c("DateTime_UTC", "AirPres_kPa"))

# Use this to convert water pressure data to water level data using sensor offset data:
# Note, as of 11/10/2019 this is a fake dataset that needs to be updated based on the field notes
sensorOffsets <- read.csv("data/raW/sensorOffsets.csv", header=T)

# List of all the variables I might want from a given site:
variables <- c("DO_mgL", "satDO_mgL","Level_m",  "WaterPres_kPa", "WaterTemp_C")

getSPdata <- function(site, start_date, end_date, variables, airP, sensorOffsets){
  site_code <- paste('NC',site,sep='_')
  dat <- request_data(site_code, start_date, end_date, variables)
  # remove bad data
    w <- which(dat$data$flagtype== "Bad Data" | dat$data$flagtype=="Questionable")
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

for(site in SPsites){
  dat <- getSPdata(site, start_date, end_date, variables, airP, sensorOffsets)
  write.csv(dat, file = paste0("data/raw/",site,".csv"),row.names = F)
  dat$site <- rep(site, nrow(dat))
  alldat <- bind_rows(alldat, dat)
}


write.csv(alldat, file = "data/raw/2019SPsites.csv", row.names = F)



#############################################################################
# Interannual comparison at NHC

start_date <- ymd_hms("2017-01-01 05:00:00 UTC") # for interannual NHC
end_date <- ymd_hms("2020-01-01 05:00:00 UTC")

ZQdat <- read_csv("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/NHC Metabolism 2019-2020/Code/siteData/NC_streampulseZQ_data.csv")

sitecode <- "NC_NHC"
variables <- c("DO_mgL", "AirPres_kPa",  "WaterPres_kPa", "WaterTemp_C")

dat<- request_data(sitecode, start_date, end_date, variables)
dat <- dat$data

w<- which(dat$flagtype=="Bad Data")
dat$value[w]<- NA

dat <- select(dat, DateTime_UTC, variable, value)%>% spread(variable, value)

dat$level_m <- calculate_level(dat, sensorOffsets$offset.m[sensorOffsets$site=="NHC"])
dat$discharge_cms <- calculate_discharge(dat,Z=ZQdat$level_m[ZQdat$site=="NHC"], Q=ZQdat$discharge_cms[ZQdat$site=="NHC"])
dat$satDO_mgL <- calc_DO_sat(dat$WaterTemp_C, dat$AirPres_kPa*10)
dat$persatDO <- dat$DO_mgL/dat$satDO_mgL

write.csv(dat, "data/raw/NHCdat.csv", row.names = F)




# Calculate daily means and rolling averages
alldat$site <- factor(alldat$site, levels = SPsites)

alldat$DateTime <- with_tz(alldat$DateTime_UTC, tzone="EST")
alldat$Date <- as.Date(alldat$DateTime, tz="EST")

DailyDO <- alldat %>% dplyr::group_by(site, Date) %>% 
  dplyr::summarise(mean.DO_mgL = mean(DO_mgL, na.rm=T),
                   min.DO_mgL = min(DO_mgL, na.rm=T), 
                   max.DO_mgL = max(DO_mgL, na.rm=T))
DailyDO <- DailyDO %>% 
  group_by(site) %>% 
  mutate(roll_mean7 = rollmean(mean.DO_mgL, 7, na.pad = T),
         roll_mean_min7 = rollmean(min.DO_mgL, 7, na.pad=T))
DailyDO <- do.call(data.frame,lapply(DailyDO, function(x) replace(x, is.infinite(x),NA)))


# Plot daily means and rolling means to compare to EPA regulation criteria
ggplot(DailyDO, aes(Date, mean.DO_mgL, group = 1)) +
  geom_ribbon(aes(ymin = min.DO_mgL, ymax = max.DO_mgL, alpha=.5))+
  #geom_line(color = "black") +
  geom_line(aes(y=min.DO_mgL), color = "black")+
  geom_line(aes(y = roll_mean7), color = "black", linetype="dashed") +
  geom_hline(yintercept = 3, color="red")+
  geom_hline(yintercept=6, color = "red", linetype="dashed")+
  facet_wrap(~site, nrow=3)+
  ylab("DO (mg/L)")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")

n.mud<- nrow(DailyDO[DailyDO$site=="Mud",])
n3.mud <- length(which(DailyDO[DailyDO$site=="Mud",]$min.DO_mgL < 3))
n6.mud<- length(which(DailyDO[DailyDO$site=="Mud",]$roll_mean7 < 6))

n.NHC<- nrow(DailyDO[DailyDO$site=="NHC",])
n3.NHC <- length(which(DailyDO[DailyDO$site=="NHC",]$min.DO_mgL < 3))
n6.NHC<- length(which(DailyDO[DailyDO$site=="NHC",]$roll_mean7 < 6))

n.UNHC<- nrow(DailyDO[DailyDO$site=="UNHC",])
n3.UNHC <- length(which(DailyDO[DailyDO$site=="UNHC",]$min.DO_mgL < 3))
n6.UNHC<- length(which(DailyDO[DailyDO$site=="UNHC",]$roll_mean7 < 6))

n3.UNHC/n.UNHC



par(mfrow = c(3,1), mar = c(0,0,0,0), oma = c(4,4,1,1))



for(i in 1:3){
  plot(DailyDO$Date[DailyDO$site==SPsites[i]],DailyDO$min.DO_mgL[DailyDO$site==SPsites[i]], 
      ylim = c(0,15),type = "l", xaxt = "n")
 # polygon(c(DailyDO$Date[DailyDO$site==SPsites[i]],rev(DailyDO$Date[DailyDO$site==SPsites[i]])),
  #        c(DailyDO$min.DO_mgL[DailyDO$site==SPsites[i]],rev(DailyDO$max.DO_mgL[DailyDO$site==SPsites[i]])),
  #        col = "red")
  lines(DailyDO$Date[DailyDO$site==SPsites[i]],DailyDO$roll_mean7[DailyDO$site==SPsites[i]], lty = 2)
  abline(h=6, col = "red", lty = 2)
  abline(h=3, col = "red")
  if(i==3) axis(1)
}
