####################
# Calculate metrics from DO timeseries
# A Carter
# 2020 Apr 16

#install.packages("streamMetabolizer", dependencies=TRUE, 
#                 repos=c("https://owi.usgs.gov/R","https://cran.rstudio.com"))

library(lubridate)
library(dplyr)
library(streamMetabolizer)
library(zoo)

# Read in compiled DO data
DOdat <- read.csv("data/raw/allSites_20180706.csv", header = T, stringsAsFactors = F)
DOdat$DateTime_UTC <- ymd_hms(DOdat$DateTime_UTC)
NHCsites2018 <- c("Mud","MC751","MC3","MC2","MC1","UNHC","NHC",
                  "NHC5","NHC4","NHC3","NHC2","NHC1")
startdate <- ymd_hms("2018-06-12 13:45:00 UTC")
enddate <- ymd_hms("2018-07-06 02:45:00 UTC")
#Convert to local time
DOdat$DateTime <- with_tz(DOdat$DateTime_UTC, tzone="EST")
DOdat$Date <- as.Date(DOdat$DateTime, tz="EST")
DOdat$Hour <- hour(DOdat$DateTime)

# Set netagive DO values to zero. Do I need to refine this?
w <-which(DOdat$DO_mgL<0)
DOdat[w, c(3,7)] <- 0
DOdat[DOdat$site=="UNHC"&DOdat$DateTime==ymd_hms("2018-06-27 11:00:00", tz="EST"), c(3,7)]<- NA
DOdat[DOdat$site=="NHC3"&DOdat$DateTime==ymd_hms("2018-06-25 13:00:00", tz="EST"), c(3,7)]<- NA
DOdat[DOdat$site=="NHC3"&DOdat$DateTime==ymd_hms("2018-06-25 12:45:00", tz="EST"), c(3,7)]<- NA


DOdat$site <- factor(DOdat$site, levels = NHCsites2018)
DOdat <- DOdat[which(DOdat$DateTime>=startdate),]
DOdat<- DOdat[which(DOdat$DateTime<=enddate),]
DOdat.mintimes <- DOdat[which(DOdat$Hour %in% c(0,1,2,3,4,5,6,7)),]
DOdat.mintimes <- DOdat.mintimes[which(DOdat.mintimes$Date >as.Date("2018-06-23")),]
DOdat.maxtimes <- DOdat[which(DOdat$Hour %in% c(11,12,13,14,15,16,17,18,19,20)),]
DOdat.maxtimes <- DOdat.maxtimes[which(DOdat.maxtimes$Date >as.Date("2018-06-23")),]

DOdatdaily<- DOdat %>% group_by(site, Date)%>%
  summarize(mean.persatDO = mean(persatDO, na.rm=T),
            sd.persatDO = sd(persatDO, na.rm=T))
DOdatdaily.mintimes <- DOdat.mintimes %>% group_by(site,Date)%>%
  summarize(min.persatDO = min(persatDO, na.rm=T), 
            time.min.persatDO = DateTime[which(persatDO==min(persatDO, na.rm=T))[ceiling(length(which(persatDO==min(persatDO, na.rm=T)))/2)]])
DOdatdaily <- full_join(DOdatdaily, DOdatdaily.mintimes, by = c("site","Date"))          
DOdatdaily.maxtimes <- DOdat.maxtimes %>% group_by(site,Date)%>%
  summarize(max.persatDO = max(persatDO, na.rm=T), 
            time.max.persatDO = DateTime[which(persatDO==max(persatDO, na.rm=T))[ceiling(length(which(persatDO==max(persatDO, na.rm=T)))/2)]])
DOdatdaily <- full_join(DOdatdaily, DOdatdaily.maxtimes, by = c("site","Date"))          

DOdatdaily$amp.persatDO <- DOdatdaily$max.persatDO-DOdatdaily$min.persatDO

par(mar = c(1,1,0,0), mfrow = c(4,3))
for(i in 1:12){
  plot(DOdat$DateTime[DOdat$site==NHCsites2018[i]], DOdat$persatDO[DOdat$site==NHCsites2018[i]],type="l", lwd="2")
  points(DOdatdaily$time.min.persatDO[DOdatdaily$site==NHCsites2018[i]],DOdatdaily$min.persatDO[DOdatdaily$site==NHCsites2018[i]],
         col = "red", pch = 19)
}




DOmetrics <- data.frame(site=NHCsites2018,
                        meanDOpersat.6.25 = 999,
                        minDOpersat.6.25 = 999,
                        ampDOpersat.6.25=999,
                        pulseDOpersat=999,
                        EOD.perday = 999,
                        dAdT.perday3=999,
                        dAdT.perday4=999,
                        dAdT.perday5=999)

stormDate <- as.Date("2018-06-26")
startDate <- as.POSIXct("2018-06-24 00:00:00")
endDate <- as.POSIXct("2018-07-04 00:00:00")

par(mar = c(1,1,0,0), mfrow = c(4,3))
for(i in 1:12){
  plot(DOdat$DateTime[DOdat$site==NHCsites2018[i]], DOdat$persatDO[DOdat$site==NHCsites2018[i]],type="l", lwd="2", 
       xlim = c(startDate,endDate), ylim = c(0,1.2), ylab = "", xlab = "")
  points(DOdatdaily$time.min.persatDO[DOdatdaily$site==NHCsites2018[i]],DOdatdaily$min.persatDO[DOdatdaily$site==NHCsites2018[i]],
         col = "brown3", pch = 19, cex=2)
  points(DOdatdaily$time.max.persatDO[DOdatdaily$site==NHCsites2018[i]],DOdatdaily$max.persatDO[DOdatdaily$site==NHCsites2018[i]],
         col = "steelblue", pch = 19, cex=2)
  w<-which(DOdat$Hour==0)
  abline(v=DOdat$DateTime[w], col="grey80")
}


for(i in 1:12){
  # daily metrics
  DO <- DOdatdaily[DOdatdaily$site==NHCsites2018[i]&DOdatdaily$Date==(stormDate-1),]
  DOmetrics$meanDOpersat.6.25[i] <- DO$mean.persatDO
  DOmetrics$minDOpersat.6.25[i]<- DO$min.persatDO
  DOmetrics$ampDOpersat.6.25[i]<- DO$amp.persatDO
  # Storm Pulse
  DO <- DOdatdaily[DOdatdaily$site==NHCsites2018[i]&DOdatdaily$Date %in% c(stormDate,stormDate-1),]
  DOmetrics$pulseDOpersat[i] <- DO$max.persatDO[2]-DO$min.persatDO[1]
  # EOD
  DO <- DOdatdaily[DOdatdaily$site==NHCsites2018[i]&DOdatdaily$Date %in% c(stormDate, stormDate+1, stormDate+2, stormDate+3, stormDate+4),]
  DO$min.persatDO[DO$Date==stormDate] <- DO$max.persatDO[DO$Date==stormDate]
  DO$time.min.persatDO[DO$Date==stormDate]<- DO$time.max.persatDO[DO$Date==stormDate]
  m<- lm(min.persatDO~time.min.persatDO, data=DO)
  DOmetrics$EOD.perday[i]<- -m$coef[2]*60*60*24 # units of DO%sat/day
  # dAmplitude/dtime
  DO <- DOdatdaily[DOdatdaily$site==NHCsites2018[i]&DOdatdaily$Date %in% c(stormDate-1, stormDate+1, stormDate+2, 
                                                                           stormDate+3),]
  DO$amp.persatDO <- DO$amp.persatDO/DO$amp.persatDO[DO$Date==(stormDate-1)]
  DO <- DO[-1,]
  DO$d <- seq(1:nrow(DO))
  m <- lm(amp.persatDO~d, data=DO)
  DOmetrics$dAdT.perday3[i] <- m$coef[2]
  
  # dAmplitude/dtime
  DO <- DOdatdaily[DOdatdaily$site==NHCsites2018[i]&DOdatdaily$Date %in% c(stormDate-1, stormDate+1, stormDate+2, 
                                                                           stormDate+3, stormDate+4),]
  DO$amp.persatDO <- DO$amp.persatDO/DO$amp.persatDO[DO$Date==(stormDate-1)]
  DO <- DO[-1,]
  DO$d <- seq(1:nrow(DO))
  m <- lm(amp.persatDO~d, data=DO)
  DOmetrics$dAdT.perday4[i] <- m$coef[2]

# dAmplitude/dtime
DO <- DOdatdaily[DOdatdaily$site==NHCsites2018[i]&DOdatdaily$Date %in% c(stormDate-1, stormDate+1, stormDate+2, 
                                                                         stormDate+3, stormDate+4, stormDate+5),]
DO$amp.persatDO <- DO$amp.persatDO/DO$amp.persatDO[DO$Date==(stormDate-1)]
DO <- DO[-1,]
DO$d <- seq(1:nrow(DO))
m <- lm(amp.persatDO~d, data=DO)
DOmetrics$dAdT.perday5[i] <- m$coef[2]
}
DOmetrics[,2:9] <- round(DOmetrics[,2:9]*100, 1)
DOmetrics

write_csv(DOmetrics, "data/DOtimeseriesMetrics_04-16-20.csv")
for(i in 1:12){
  # dAmplitude/dtime
  DO <- DOdatdaily[DOdatdaily$site==NHCsites2018[i]&DOdatdaily$Date %in% c(stormDate-1, stormDate+1, stormDate+2, 
                                                                           stormDate+3, stormDate+4, stormDate+5),]
  DO$amp.persatDO <- DO$amp.persatDO/DO$amp.persatDO[DO$Date==(stormDate-1)]
  DO <- DO[-1,]
  DO$d <- seq(1:nrow(DO))
  m <- lm(amp.persatDO~d, data=DO)
  plot(DO$d, DO$amp.persatDO, col="brown3", pch=19, cex=1.5) 
  abline(m$coef[1], m$coef[2])
  
  DO <- DOdatdaily[DOdatdaily$site==NHCsites2018[i]&DOdatdaily$Date %in% c(stormDate-1, stormDate+1, stormDate+2, 
                                                                           stormDate+3, stormDate+4),]
  DO$amp.persatDO <- DO$amp.persatDO/DO$amp.persatDO[DO$Date==(stormDate-1)]
  DO <- DO[-1,]
  DO$d <- seq(1:nrow(DO))
  m <- lm(amp.persatDO~d, data=DO)

  abline(m$coef[1], m$coef[2], lty=2)
  
  DO <- DOdatdaily[DOdatdaily$site==NHCsites2018[i]&DOdatdaily$Date %in% c(stormDate-1, stormDate+1, stormDate+2, 
                                                                           stormDate+3),]
  DO$amp.persatDO <- DO$amp.persatDO/DO$amp.persatDO[DO$Date==(stormDate-1)]
  DO <- DO[-1,]
  DO$d <- seq(1:nrow(DO))
  m <- lm(amp.persatDO~d, data=DO) 
  abline(m$coef[1], m$coef[2],lty=3)
}


#######################################
#load annual data for flashiness and sagginess

dat <- read.csv("data/raw/2019SPsites.csv", header = T, stringsAsFactors = F)
dat$DateTime_UTC <- ymd_hms(dat$DateTime_UTC)
dat$DateTime <- with_tz(dat$DateTime_UTC, tz="EST")
dat$Date <- as.Date(dat$DateTime, tz = "EST")
dat <- dat[dat$site!="NHC",]


datNHC<- read.csv("data/raw/NHCdat.csv")
datNHC$DateTime_UTC<- ymd_hms(datNHC$DateTime_UTC)
datNHC$DateTime <- with_tz(datNHC$DateTime_UTC, tz="EST")
datNHC$Date <- as.Date(datNHC$DateTime, tz="EST")
datNHC$site <- "NHC"
# Summarize by day
datdaily <- dat %>% dplyr::group_by(site, Date) %>% 
  dplyr::summarise(mean.persatDO = mean(persatDO, na.rm=T),
                   min.persatDO=min(persatDO, na.rm=T),
                   max.persatDO=max(persatDO,na.rm=T),
                   sd.persatDO = sd(persatDO, na.rm=T)) 
datdaily$mean.persatDO <- na.approx(datdaily$mean.persatDO)

datdailyNHC <- datNHC %>% dplyr::group_by(site, Date) %>% 
  dplyr::summarise(mean.persatDO = mean(persatDO, na.rm=T),
                   min.persatDO=min(persatDO, na.rm=T),
                   max.persatDO=max(persatDO,na.rm=T),
                   sd.persatDO = sd(persatDO, na.rm=T)) 
datdailyNHC$mean.persatDO <- na.approx(datdailyNHC$mean.persatDO)

datdaily <- rbind(datdaily, datdailyNHC)
datdaily$year <- year(datdaily$Date)

annualSites <- c("Mud","MC751", "UNHC", "NHC")
# Calculate RBI for each time series
# m is a DO timeseries, t is the corresponding time intervals.
# The resulting metric is in units of DO %sat/time
# Make sure that time intervals are 1 hour for consistency
RBIcalc <- function(m,t){
  l <- abs(diff(m, na.rm=T))
  RBI <- sum(l, na.rm=T)/as.numeric((t[length(t)]-t[1]))
  RBI
}

sag_calc <- function(DO, t){
  DO <- na.approx(DO)
  DO <- sort(DO, decreasing=T)
  index <- seq(1:length(DO))
  freq <- index/(length(DO))
  l <- diff(freq)[1]
  DO.m <- rep(1,length(DO)-1)
  for(i in 1:(length(DO)-1)){
    DO.m[i] <- (DO[i]+DO[i+1])/2
  }
  area <- sum(DO.m)*l
  area
}
DOannualmetrics <- datdaily %>% dplyr::group_by(site,year) %>%
  dplyr::summarise(meanDO = mean(mean.persatDO, na.rm=T),
                   minDO = min(mean.persatDO, na.rm=T),
                   maxDO = max(mean.persatDO, na.rm=T))
DOannualmetrics$flashiness <- 999
DOannualmetrics$sagginess <- 999

for(i in 1:nrow(DOannualmetrics)){
  site <- DOannualmetrics$site[i]
  year <- DOannualmetrics$year[i]
  DO <- datdaily[datdaily$site==site&datdaily$year==year,]
  DOannualmetrics$flashiness[i] <- RBIcalc(DO$mean.persatDO, DO$Date)
  DOannualmetrics$sagginess[i] <- sag_calc(DO$mean.persatDO, DO$Date)
}


DOannualmetrics[,3:7] <- round(DOannualmetrics[,3:7]*100,1)
write_csv(DOannualmetrics, "data/DOAnnualMetrics_04-16-20.csv")
