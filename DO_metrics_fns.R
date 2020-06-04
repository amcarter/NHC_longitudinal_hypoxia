#####################
# Functions to calculate DO metrics
#   Storm identification - coming soon
#   Storm events:
#     slopes of minima
#     slopes of amplitudes
#   Annual:
#     DO undersaturation curves
#     Day:night hypox ratio

# A Carter 2020 05 22

library(lubridate)
library(dplyr)
library(readr)
library(streamMetabolizer)

dat <- read_csv("data/raw/NHCdat.csv")
sites <- read.csv("NHC_map/NC_synopticSamplingSites.csv", header=T, stringsAsFactors = F)

stormdate <- as.Date("2018-06-26")
DO <- dat$DO_mgL
dt <- dat$DateTime_UTC
local_dt <- with_tz(dt, tz="EST")


###########################################################################
# Slope of mins and maxes and amplitude recovery after a storm

calc_DO_storm_recovery<- function(local_dt,DO, stormdate,minwin=0:6, maxwin=9:19, rec_days=5){
  dat<- data.frame(dt = local_dt,                  
                   DO=DO,                    
                   date=as.Date(floor_date(local_dt, "days")),
                   hour=hour(local_dt))
  xlims <- c(stormdate-2, stormdate+6)
  tz <- attr(local_dt, "tzone")
  xlims <- as.POSIXct(paste(as.character(xlims),"00:00:00",sep=" "), tz=tz)
  dat<- dat[dat$dt<xlims[2]&dat$dt>=xlims[1],]
  
  
  #Calculate post storm peak
  peak_DO <- max(dat$DO[dat$date == stormdate], na.rm=T)
  peak_time <- dat[dat$date == stormdate,]$dt[
    which(dat$DO[dat$date == stormdate]==peak_DO)][1]
 
  # data frame of mins and maxes:
  datmin <- dat[dat$hour %in% minwin,]
  datmax <- dat[dat$hour %in% maxwin,]
 
  daily<- data.frame(date=stormdate, min_time=peak_time, min_DO=peak_DO, max_time=peak_time, max_DO=peak_DO)
  days <- unique(dat$date)
  for(i in 1:length(days)){
    tmp <-  datmin[datmin$date==days[i],]
    min_DO <- min(tmp$DO, na.rm=T)
    min_time <- tmp$dt[which(tmp$DO==min_DO)][ceiling(length(which(tmp$DO==min_DO))/2)] # find the middle time point of mins
  
    tmp <-  datmax[datmax$date==days[i],]
    max_DO <- max(tmp$DO, na.rm=T)
    max_time <- tmp$dt[which(tmp$DO==max_DO)][ceiling(length(which(tmp$DO==max_DO))/2)] 
  
    newrow <- data.frame(date=days[i], min_time=min_time, min_DO=min_DO, 
                         max_time=max_time, max_DO=max_DO)
    daily<- bind_rows(daily, newrow)
  }
 
  daily$amp <- daily$max_DO-daily$min_DO
  
  # remove mins and maxes before the peak
  preamp <- daily$amp[daily$date==stormdate-1]
  daily<- daily[-which(daily$max_time<peak_time),]
  w<- which(daily$min_time<peak_time)
  if(length(w)!=0){
    daily<- daily[-w,]
  }
  daily<- daily[1:(rec_days+1),]
  # if(nrow(daily)<6){
  #   stop("not enough days to calculate metric, need 5 days post storm")}
  
  minm<- lm(min_DO~min_time, data=daily)
  maxm<- lm(max_DO~max_time, data=daily)
  min_slope<- minm$coefficients[2]
  min_int<- minm$coefficients[1]
  max_slope<- maxm$coefficients[2]
  max_int<- maxm$coefficients[1]
  min_r2<-summary(minm)$adj
  max_r2<-summary(maxm)$adj
  dDO.day_min <- min_slope*60*60*24 # convert to dDO/day
  dDO.day_max <- max_slope*60*60*24 # convert to dDO/day
  
  plot(dat$dt, dat$DO, col="grey40", type="l",lwd=1.2,
       xaxt="n",yaxt="n",xlab="",ylab="",ylim = c(0,130))
  abline(min_int,min_slope, lwd=2, lty=2, col="grey50")
  abline(min_int,min_slope, lwd=2, lty=2, col=alpha("brown3",.6))
  abline(max_int,max_slope, lwd=2, lty=2, col="grey50")
  abline(max_int,max_slope, lwd=2, lty=2, col=alpha("steelblue", alpha=.6))
  points(daily$min_time, daily$min_DO, col="brown3", pch=19)
  points(daily$max_time, daily$max_DO, col="steelblue", pch=19)
  points(daily$max_time[1], daily$max_DO[1], pch=19, cex=1.2)
 
  # mtext(paste0("min r2 = ",round(min_r2,2),
  #              "    max r2 = ",round(max_r2,2)),1, -2)
  params <- data.frame(dDO.day_min=dDO.day_min,
                       dDO.day_max=dDO.day_max,
                       pre_amplitude=preamp,
                       amp_recovery_percent.day=-(dDO.day_min-dDO.day_max),
                       r2.adj_min=min_r2,
                       r2.adj_max=max_r2)
  
  DOfit <- list(dat=dat, daily=daily,params=params)
  return(DOfit)
}



calc_night_hypoxia<- function(UTC_dt, DO, threshold=.5, lat, long){
  dat<- data.frame(dt = UTC_dt,
                   DO = DO)
  dat$solar_time <- convert_UTC_to_solartime(UTC_dt, long, "apparent solar")
  dat$insolation <- calc_solar_insolation(dat$solar_time, lat)
  dat$light<- 1
  dat$light[dat$insolation==0]<- 0
  dat$hypox<-0
  dat$hypox[dat$DO<=threshold]<-1
  dat<- dat[!is.na(dat$DO),]
  
  ts_length <-nrow(dat)/96 
  per_hypox <- sum(dat$hypox)/nrow(dat)
  per_night_hypox <- sum(dat$hypox[dat$light==0])/nrow(dat[dat$light==0,])
  per_day_hypox <- sum(dat$hypox[dat$light==1])/nrow(dat[dat$light==1,])
  
  tmp<- dat[dat$hypox==1,]
  fraction_night_hypoxia <- 1-sum(tmp$light)/nrow(tmp)
  
  prob_night_given_hypoxia <- per_night_hypox/per_day_hypox
  
  # plot(dat$solar_time, dat$DO, type="l", ylim = c(0,1.2))
  # abline(v=dat$dt[dat$light==0], col=alpha("black",.01))
  # abline(h=threshold, lty=2, col="red")
  out <- data.frame(ts_days=ts_length, 
                    per_hypox=per_hypox, 
                    per_night_hypox=per_night_hypox,
                    per_day_hypox=per_day_hypox,
                    prob_night_given_hypoxia=prob_night_given_hypoxia, 
                    hypoxia_threshold=threshold)
  return(out)
  }
