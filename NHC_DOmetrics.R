####################
# Summary statistics for NHC 2017-2019
# A Carter
# 2020 Apr 24

library(lubridate)
library(dplyr)
library(readr)
library(zoo)

# Read in compiled DO data


dat<- read.csv("data/raw/NHCdat.csv")
dat$DateTime_UTC<- ymd_hms(dat$DateTime_UTC)
dat$DateTime <- with_tz(dat$DateTime_UTC, tz="EST")
dat$month <- month(dat$DateTime)
dat$Date <- as.Date(dat$DateTime, tz="EST")
dat$year <- year(dat$Date)
# measurements below 50% saturation
w <- which(dat$persatDO<.5)
length(w)/length(!is.na(dat$persatDO))

datlowDO <- dat[w,]
# percent of low DO that occurs in October
length(which(datlowDO$month==10))/nrow(datlowDO)
dat2019 <- dat[dat$year==2019,]
datdaily <-  dat2019 %>% dplyr::group_by(Date) %>% 
  dplyr::summarise(mean.persatDO = mean(persatDO, na.rm=T),
                   min.persatDO=min(persatDO, na.rm=T),
                   max.persatDO=max(persatDO,na.rm=T),
                   sd.persatDO = sd(persatDO, na.rm=T), 
                   amp.persatDO = max.persatDO-min.persatDO,
                   Q = mean(discharge_cms, na.rm=T)) 
datdaily$amp.persatDO[datdaily$amp.persatDO==-Inf] <- NA
datdaily$Q <- na.approx(datdaily$Q, na.rm=F)
# datdaily$mean.persatDO <- na.approx(datdaily$mean.persatDO)

# Filter out all of the storm dates

plot(dat2019$DateTime[n:m], dat2019$discharge_cms[n:m], ylim = c(.1,100),log = "y", type =  "l")
abline(h=quantile(dat2019$discharge_cms, .9, na.rm=T))

tmp<- diff(datdaily$Q)
tmp[tmp<0]<- 0
plot(datdaily$Date[301:362], tmp[300:361], type = "l", ylim = c(0,2*q.9))
abline(h=q.9)

stormdates <- as.Date(c("2019-01-13","2019-01-20","2019-01-24","2019-02-13",
                        "2019-02-16","2019-02-18","2019-02-21","2019-02-24",
                        "2019-03-21","2019-04-06","2019-04-09","2019-04-13",
                        "2019-04-20","2019-06-08","2019-06-19","2019-07-23",
                        "2019-08-05","2019-08-14","2019-10-20","2019-11-24",
                        "2019-12-01","2019-12-14", "2019-12-24"))

NHCmetrics <- data.frame(
  stormDate = stormdates,
  stormAmp=999,
  EOD = 999,
  dAdT = 999)

for(i in 1:length(stormdates)){
  stormDate<- stormdates[i]
  # Storm Pulse
  DO <- datdaily[datdaily$Date %in% c(stormDate,stormDate-1),]
  NHCmetrics$stormAmp[i] <- DO$max.persatDO[2]-DO$min.persatDO[1]
  # EOD
  DO <- datdaily[datdaily$Date %in% c(stormDate, stormDate+1, stormDate+2, stormDate+3, stormDate+4),]
  DO$min.persatDO[DO$Date==stormDate] <- DO$max.persatDO[DO$Date==stormDate]
  DO$time.min.persatDO[DO$Date==stormDate]<- DO$time.max.persatDO[DO$Date==stormDate]
  m<- lm(min.persatDO~time.min.persatDO, data=DO)
  DOmetrics$EOD.perday[i]<- -m$coef[2]*60*60*24 # units of DO%sat/day
  # dAmplitude/dtime
  DO <- DOdatdaily[DOdatdaily$site==NHCsites2018[i]&DOdatdaily$Date %in% c(stormDate-1, stormDate+1, stormDate+2, 
                                                                           stormDate+3, stormDate+4),]
  DO$amp.persatDO <- DO$amp.persatDO/DO$amp.persatDO[DO$Date==(stormDate-1)]
  DO <- DO[-1,]
  DO$d <- seq(1:nrow(DO))
  m <- lm(amp.persatDO~d, data=DO)
  DOmetrics$dAdT.perday[i] <- m$coef[2]
}