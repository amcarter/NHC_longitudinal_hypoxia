####################
# Calculate metrics from DO timeseries
# A Carter
# 2019 Nov 10

#install.packages("streamMetabolizer", dependencies=TRUE, 
#                 repos=c("https://owi.usgs.gov/R","https://cran.rstudio.com"))

library(lubridate)
library(dplyr)
library(streamMetabolizer)
library(zoo)
library(scales)

setwd(hypox_projdir)
source("DO_metrics_fns.R")
# Read in compiled DO data

sites <- read.csv("NHC_map/NC_synopticSamplingSites.csv", header=T, stringsAsFactors = F)

DOdat <- read.csv("data/raw/allSites_20180706.csv", header = T, stringsAsFactors = F)
DOdat$DateTime_UTC <- ymd_hms(DOdat$DateTime_UTC)
NHCsites2018 <- c("Mtrib","MC751","MC3","MC2","MC1","UNHC","NHC",
                  "NHC5","NHC4","NHC3","NHC2","NHC1")
DOdat$site[DOdat$site=="Mud"]<- "Mtrib"
startdate <- ymd_hms("2018-06-12 13:45:00 UTC")
enddate <- ymd_hms("2018-07-06 02:45:00 UTC")

#DOdat <- DOdat[which(DOdat$DateTime_UTC>=startdate),]
#DOdat<- DOdat[which(DOdat$DateTime_UTC<=enddate),]

#Convert to local time
DOdat$DateTime <- with_tz(DOdat$DateTime_UTC, tzone="EST")
DOdat$Date <- as.Date(DOdat$DateTime, tz="EST")

# Set netagive DO values to zero. Do I need to refine this?
w <-which(DOdat$DO_mgL<0)
DOdat[w, c(3,7)] <- 0

# Calculate metrics:
# Annual: mean, range, night hypox ratio, num hypoxic events, max hypoxia duration, modal undersaturation, ____undersaturation
# seasonal: night hypox ratio, pr(hypox)pr(anox)
# Event: Storm pulse, oxygen demand, amplitude recovery
# daily: mean, min, amplitude

#hypoxia thresholds:
th1=.25
th2=.5

# daily metrics:
# calculated for pre-storm day, 6/25/2018

prestorm <- as.Date("2018-06-25")
DO_day <- DOdat[DOdat$Date==prestorm,] %>%group_by(site) %>%
  summarise(mean_DO_persat=mean(persatDO, na.rm=T),
            min_DO_persat=min(persatDO, na.rm=T),
            max_DO_persat=max(persatDO, na.rm=T))

DO_day$amp_DO_persat <- DO_day$max_DO_persat-DO_day$min_DO_persat
colnames(DO_day)<- c("site",  paste0(colnames(DO_day[,2:5]),".6.25"))

# multiple days/event:
# calculated from whole considered timeperiod
stormdate<- as.Date("2018-06-26")
DO_event<- data.frame()
for(site in NHCsites2018){
  dat<- DOdat[DOdat$site==site,]
  # probability of hypoxia
  cdf<- ecdf(dat$persatDO)
  pr_anox <- cdf(0)
  pr_hypox_th1<- cdf(th1)
  pr_hypox_th2<- cdf(th2)
  # ratio of night hypox/total hypox
  night_hypox_ratio <- calc_night_hypoxia(dat$DateTime, dat$persatDO, threshold = .5, 
                                          lat=sites$Lat[sites$site==site], 
                                          long=sites$Long[sites$site==site])$night_hypoxia_ratio
  max_hypox_interval_days <- hypoxia_durations(dat$persatDO)$max_hrs/24
  max_anox_interval_days <- hypoxia_durations(dat$persatDO, threshold = 0)$max_hrs/24
  # storm event metrics
  storm <- calc_DO_storm_recovery(dat$DateTime, dat$persatDO, stormdate)$params
  
  tmp <- data.frame(site=site, days=nrow(dat)/96,
                    pr_anox=pr_anox, pr_hypox_th1=pr_hypox_th1, pr_hypox_th2=pr_hypox_th2, 
                    night_hypox_ratio=night_hypox_ratio,
                    max_anox_interval_days=max_anox_interval_days,
                    max_hypox_interval_days=max_hypox_interval_days,
                    dDO.day_min=storm$dDO.day_min, amp_recovery.day=storm$amp_recovery_percent.day, stormpulse=storm$stormpulse)
  DO_event <- rbind(DO_event, tmp)
}

sites$dist_upstream <- sites$pathlength-sites$pathlength[sites$site=="NHC1"]
DO_metrics <- full_join(DO_day, DO_event, by="site")
DO_metrics <- left_join(DO_metrics, sites[,c(1,9,10)])



write_csv(DO_metrics, "data/DO_metrics_2020_06_03.csv")

###########################################################################
#calculate annual metrics

#######################################
#load annual data for flashiness and sagginess

annual_dat <- read.csv("data/raw/2019SPsites.csv", header=T, stringsAsFactors = F)
annual_dat$DateTime_UTC <- ymd_hms(annual_dat$DateTime_UTC)
annual_dat$DateTime <- with_tz(annual_dat$DateTime_UTC, tz="EST")
annual_dat$date <- as.Date(annual_dat$DateTime, tz = "EST")
annual_dat$year <- year(annual_dat$date)

yearsites <- c("NHC","UNHC","Mud","MC751")
sites$site[sites$site=="Mtrib"]<- "Mud"
DO_annual<- data.frame()
for(site in yearsites){
  lat = sites$Lat[sites$site==site]
  long= sites$Long[sites$site==site]
    
  siteyears=c(2017,2018,2019)
  if(site=="MC751"){siteyears=2019}
  
  for(i in siteyears){
    dat <- annual_dat[annual_dat$site==site&annual_dat$year==i,]
    meanDO <- mean(dat$persatDO, na.rm=T)
    minDO <- min(dat$persatDO, na.rm=T)
    maxDO <- max(dat$persatDO, na.rm=T)
    night <- calc_night_hypoxia(dat$DateTime, dat$persatDO, .5, lat, long)
    intervals <- hypoxia_durations(dat$persatDO)
    
    ### duration statistics
    # Summarize by day
    datdaily <- dat %>% dplyr::group_by(date) %>% 
      dplyr::summarise(mean.persatDO = mean(persatDO, na.rm=T),
                     min.persatDO=min(persatDO, na.rm=T),
                     max.persatDO=max(persatDO,na.rm=T))
    datdaily$mean.persatDO[is.infinite(datdaily$mean.persatDO)]<- NA
    modeDO_daily <- calcmode(round(datdaily$mean.persatDO,2))
    meanDO_daily <- mean(datdaily$mean.persatDO, na.rm=T)
    varDO_daily<- var(datdaily$mean.persatDO, na.rm=T)
    mu <- log(meanDO_daily^2/sqrt(meanDO_daily^2+varDO_daily))
    sigma <- sqrt(log(1+varDO_daily/meanDO_daily^2))
    
    hist(datdaily$mean.persatDO, main="", xlab="", xaxt="n", yaxt="n", ylab="n")
    par(new=T)
    
    x <- seq(0,1, by=.05)
    y <- dlnorm(x, meanlog=(meanDO_daily), sdlog=.4)
    plot(x,y)
    plot(ecdf(datdaily$mean.persatDO), main=paste(site,i))
    
    #datdaily$mean.persatDO <- na.approx(datdaily$mean.persatDO)

    
    tmp <- data.frame(site=site, year=i, min_DO_persat=minDO,max_DO_persat=maxDO, mean_DO_persat=meanDO,
                      mode_DO_daily_persat=modeDO_daily,
                      night_hypox_ratio=night$night_hypoxia_ratio, per_hypox=night$per_hypox,
                      max_hypox_interval_days=intervals$max_hrs/24, num_hypox_intervals = length(intervals$interval_lengths)
                      )
    DO_annual <- rbind(DO_annual, tmp)
  }
}
  
  

write_csv(DO_annual, "data/DO_annual_metrics_2020_06_03.csv")



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























# YSIdat <- read_csv("data/FieldConcFluxes_NHCMC_20180702.csv")%>%
#   select(site, SOM = SedOrgM.per, ysi.DO_mgL = DO.mgL.ysi, ysi.persatDO = DO.persat.ysi)
# YSIdat <- YSIdat[-4,]
# YSIdat$site[3] <- "MC3"
# # Create a dataframe with only the 9-5 observations
# 
# w <- c(which(hour(DOdat$DateTime) %in% 9:16),
#   which(hms(strftime(DOdat$DateTime,format = "%H:%M:%S"))==hms("05:00:00")))
# DOdatWD <- DOdat[w,]

# Calculate statistics for percent of time hypoxic
# Set Thresholds:
th1 <- .50
th2 <- .25

hypoDOdat <- DOdat
hypoDOdat$hypox.th2 <- rep(NA, nrow(DOdat))
hypoDOdat$hypox.th2[which(DOdat$DO_mgL>th2)]<-0
hypoDOdat$hypox.th2[which(DOdat$DO_mgL<=th2)]<-1

hypoDOdat$hypox.th1 <- rep(NA, nrow(DOdat))
hypoDOdat$hypox.th1[which(DOdat$persatDO>th1)]<-0
hypoDOdat$hypox.th1[which(DOdat$persatDO<=th1)]<-1

hypoDOdat$hypox0 <- rep(NA, nrow(DOdat))
hypoDOdat$hypox0[which(DOdat$persatDO>0)]<-0
hypoDOdat$hypox0[which(DOdat$persatDO<=0)]<-1

# Calculate RBI for each time series
#RBIcalc <- function(m){
#  l <- abs(diff(m, na.rm=T))
#  RBI <- sum(l, na.rm=T)/sum(m, na.rm=T)
#  RBI
#}

RBIcalc <- function(m,t){
  l <- abs(diff(m, na.rm=T))

  RBI <- sum(l, na.rm=T)/as.numeric((t[length(t)]-t[1]))
  RBI
}

# Create dataframes for 1, 7, full ts metrics:
sampledate <- ymd("2018-07-02")
DOdat.day <-DOdat[DOdat$Date==sampledate,]
DOdat.7day <- DOdat[DOdat$Date>(sampledate-days(7))&DOdat$Date<=sampledate,]
DOdat <- DOdat[DOdat$DateTime>=ymd_hms("2018-06-13 00:00:00", tz="EST")&
                 DOdat$DateTime<ymd_hms("2018-07-06 00:00:00", tz="EST"),]
#tmp <- tapply(DOdat$DO_mgL, DOdat$site, mean, na.rm=T)

DOsummaryStats <- data.frame(site = YSIdat$site,
                             ysi.persatDO = YSIdat$ysi.persatDO,
                             SOM.percent = YSIdat$SOM,
                             # calculate means and standard deviations
                             meanday.DO_mgL = tapply(DOdat.day$DO_mgL, DOdat.day$site, mean, na.rm=T),
                             mean7day.DO_mgL = tapply(DOdat.7day$DO_mgL, DOdat.7day$site, mean, na.rm=T),
                             meanall.DO_mgL = tapply(DOdat$DO_mgL, DOdat$site, mean, na.rm = T),
                             sdday.DO_mgL = tapply(DOdat.day$DO_mgL, DOdat.day$site, sd, na.rm=T),
                             sd7day.DO_mgL = tapply(DOdat.7day$DO_mgL, DOdat.7day$site, sd, na.rm=T),
                             sdall.DO_mgL = tapply(DOdat$DO_mgL, DOdat$site, sd, na.rm=T),
                             # RBI
                             RBIday.DO_mgL = tapply(DOdat.day$DO_mgL, DOdat.day$site, RBIcalc, t=DOdat.day$DateTime),
                             RBI7day.DO_mgL = tapply(DOdat.7day$DO_mgL, DOdat.7day$site, RBIcalc, t=DOdat.7day$DateTime),
                             RBIall.DO_mgL = tapply(DOdat$DO_mgL, DOdat$site, RBIcalc, t=DOdat$DateTime),
                             # # workday values:
                             # WDmean.DO_mgL = tapply(DOdatWD$DO_mgL, DOdatWD$site, mean, na.rm=T),
                             # WDmedian.DO_mgL = tapply(DOdatWD$DO_mgL, DOdatWD$site, median, na.rm=T),
                             # WDsd.DO_mgL = tapply(DOdatWD$DO_mgL, DOdatWD$site, sd, na.rm=T),
                             # WDmean.DO_persat = tapply(DOdatWD$persatDO, DOdatWD$site, mean, na.rm=T),
                             # WDmedian.DO_persat = tapply(DOdatWD$persatDO, DOdatWD$site, median, na.rm=T),
                             # WDsd.DO_persat = tapply(DOdatWD$persatDO, DOdatWD$site, sd, na.rm=T),
                             #Percent hypoxic
                             prHypox0 = tapply(hypoDOdat$hypox0, hypoDOdat$site, mean, na.rm=T),
                             prHypox50 = tapply(hypoDOdat$hypox.th1, hypoDOdat$site, mean, na.rm=T)
                             #prHypox.th2 = tapply(hypoDOdat$hypox.th2, hypoDOdat$site, mean, na.rm=T)
)
# Calculate daily means and daily work day means based around the EPA criteria for Freshwaters

DailyDO <- DOdat %>% dplyr::group_by(site, Date) %>% 
  dplyr::summarise(mean.DO_mgL = mean(DO_mgL, na.rm=T),
                   mean.DO_mgL.sd = sd(DO_mgL, na.rm=T),
                   min.DO_mgL = min(DO_mgL, na.rm=T), 
                   max.DO_mgL = max(DO_mgL, na.rm=T))
DailyDO$amp.DO_mgL <- DailyDO$max.DO_mgL- DailyDO$min.DO_mgL

# WDDailyDO <- DOdatWD %>% group_by(site, Date) %>% 
#   dplyr::summarise(WDmean.DO_mgL = mean(DO_mgL, na.rm=T),
#             WDmin.DO_mgL = min(DO_mgL, na.rm=T), 
#             WDmax.DO_mgL = max(DO_mgL, na.rm=T),
#             WDmed.DO_mgL = median(DO_mgL, na.rm=T))
# DailyDO <- cbind(DailyDO, select(WDDailyDO,-Date))%>%select(-site1)

DailyDO <- do.call(data.frame,lapply(DailyDO, function(x) replace(x, is.infinite(x),NA)))

tmp <- data.frame(ampday.DO_mgL = DailyDO[DailyDO$Date==sampledate,]$amp.DO_mgL,
                  amp7day.DO_mgL = tapply(DailyDO[DailyDO$Date<=sampledate&DailyDO$Date>(sampledate-days(7)),]$amp.DO_mgL,
                                           DailyDO[DailyDO$Date<=sampledate&DailyDO$Date>(sampledate-days(7)),]$site,
                                           mean),
                  sd.amp7day.DO_mgL = tapply(DailyDO[DailyDO$Date<=sampledate&DailyDO$Date>(sampledate-days(7)),]$amp.DO_mgL,
                                             DailyDO[DailyDO$Date<=sampledate&DailyDO$Date>(sampledate-days(7)),]$site,
                                             sd),
                  ampall.DO_mgL = tapply(DailyDO$amp.DO_mgL, DailyDO$site, mean, na.rm=T),
                  sd.ampall.DO_mgL= tapply(DailyDO$amp.DO_mgL, DailyDO$site, sd, na.rm=T),
                  minday.DO_mgL = DailyDO[DailyDO$Date==sampledate,]$min.DO_mgL,
                  min7day.DO_mgL = tapply(DailyDO[DailyDO$Date<=sampledate&DailyDO$Date>(sampledate-days(7)),]$min.DO_mgL,
                                          DailyDO[DailyDO$Date<=sampledate&DailyDO$Date>(sampledate-days(7)),]$site,
                                          mean),
                  sd.min7day.DO_mgL = tapply(DailyDO[DailyDO$Date<=sampledate&DailyDO$Date>(sampledate-days(7)),]$min.DO_mgL,
                                          DailyDO[DailyDO$Date<=sampledate&DailyDO$Date>(sampledate-days(7)),]$site,
                                          sd),
                  minall.DO_mgL = tapply(DailyDO$min.DO_mgL, DailyDO$site, mean, na.rm=T),
                  sd.minall.DO_mgL = tapply(DailyDO$min.DO_mgL, DailyDO$site, sd, na.rm=T))
                  
DOmetrics <- cbind(DOsummaryStats, tmp)                  

### This is calculating rolling means by letting data from different sites overlap
# This is a problem and should be fixed,***But*** it doesn't affect the
# date I am interested in, 7/2/18 so I am going to deal with this later
# DailyDO <- DailyDO %>% 
#   group_by(site) %>% 
#   mutate(roll_mean7 = rollmean(mean.DO_mgL, 7, na.pad = T),
#          roll_mean_min7 = rollmean(min.DO_mgL, 7, na.pad=T),
#          roll_mean7.sd = rollapply(mean.DO_mgL, 7, sd, na.pad = T))
# 
# DO20180702 <- DailyDO[DailyDO$Date==ymd("2018-07-02"),]%>% 
#   select(site,mean.DO_mgL.0702 = mean.DO_mgL, 
#          mean.DO_mgL.sd.0702=mean.DO_mgL.sd,
#          min.DO_mgL.0702=min.DO_mgL,
#          roll_mean7.0702=roll_mean7,
#          roll_mean7.sd.0702 = roll_mean7.sd)
# DOmetrics <- left_join(DOsummaryStats, DO20180702, by = "site")



write.csv(DOmetrics, file = "data/DOtimeseriesSummary.csv", row.names = F)

write.csv(DailyDO, file = "data/DailyDOsummary.csv", row.names = F)

# Plot daily means and rolling means to compare to EPA regulation criteria
library(ggplot2)
ggplot(DailyDO, aes(Date, mean.DO_mgL, group = 1)) +
  geom_ribbon(aes(ymin = min.DO_mgL, ymax = max.DO_mgL, alpha=.05))+
  #geom_line(color = "black") +
  geom_line(aes(y=min.DO_mgL), color = "black")+
  geom_line(aes(y = roll_mean7), color = "black", linetype="dashed") +
  geom_hline(yintercept = 3, color="red")+
  geom_hline(yintercept=6, color = "red", linetype="dashed")+
  facet_wrap(~site)+
  ylab("DO (mg/L)")+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")
  

par(mfrow = c(4, 3), mar = c(0,0,0,0))
for(i in length(NHCsites2018)){
  plot(DOdat[DOdat$site==NHCsites2018[i],]$DateTime_UTC,
       DOdat[DOdat$site==NHCsites2018[i],]$DO_mgL,
       xaxt=F, yaxt = F, xlab='',ylab='')
  
}








n <- length(NHCsites2018)
DOSiteSummary <- data.frame(site = NHCsites2018,
                            meanDO_mgL = as.double(rep(NA, n)),
                            sdDO_mgL = as.double(rep(NA,n)),
                            PrAnox = as.double(rep(NA, n)), 
                            PrHypox2 = as.double(rep(NA,n)),
                            PrHypox5 = as.double(rep(NA, n)),
                            dayMean = as.double(rep(NA,n)),
                            cross0 = as.double(rep(NA, n)),
                            cross2 = as.double(rep(NA, n)),
                            cross5 = as.double(rep(NA, n)),
                            amp.mean =as.double(rep(NA,n)),
                            amp.sd = as.double(rep(NA,n)),
                            RBI = as.double(rep(NA,n))) 
                            #ERnttm.mean=as.double(rep(NA,n)),
                            #ERnttm.sd=as.double(rep(NA,n)))



for (i in 1:9){
  site <- NHCsites2018[i]
  siteDO <- DOdat[[paste(site, "DO_mgL", sep=".")]]
  DOSiteSummary$meanDO_mgL[i] <- mean(siteDO,na.rm = T)
  DOSiteSummary$sdDO_mgL[i] <- sd(siteDO,na.rm=T)
  
  DOSiteSummary$RBI[i] <- RBIcalc(siteDO)  
  
  # Calculate the percent of measurements below hypoxia thresholds 
  n <- which(is.na(siteDO)==FALSE)
  n0 <- (which(siteDO<=0))
  n2 <- (which(siteDO<=2))
  n5 <- (which(siteDO<=5))
  DOSiteSummary$PrAnox[i] <- length(n0)/length(n)
  DOSiteSummary$PrHypox2[i] <- length(n2)/length(n)
  DOSiteSummary$PrHypox5[i] <- length(n5)/length(n)
  
  # Calculate the number of times each hypoxia threshold is crossed
  b0 <- siteDO
  b0[n]<- 1
  b0[n0] <- 0  
  DOSiteSummary$cross0[i]<- 1/(length(which(diff(b0)==1))/(length(n)/96)) # Calculate number of crosses per day
  b2 <- siteDO
  b2[n]<- 1
  b2[n2] <- 0  
  DOSiteSummary$cross2[i] <- 1/(length(which(diff(b2)==1))/(length(n)/96))
  b5 <- siteDO
  b5[n]<- 1
  b5[n5] <- 0  
  DOSiteSummary$cross5[i] <- 1/(length(which(diff(b5)==1))/(length(n)/96))
}

# Calculate average daily amplitude
DOdat$date <- date(DOdat$DateTime_UTC)
DOdat$hour <- lubridate::hour(DOdat$DateTime_UTC)

for(i in 1:9){
  site <- NHCsites2018[i]
  daily <- data.frame(DO=DOdat[[paste(site,"DO_mgL", sep=".")]],
                      date=DOdat$date)%>%
    group_by(date) %>%
    summarize(mean =mean(DO,na.rm = T),
                min = min(DO, na.rm = T),
                max = max(DO,na.rm = T),
                amp = (max-min))
  DOSiteSummary$amp.mean[i] <- mean(daily$amp, na.rm=T)
  DOSiteSummary$amp.sd[i]<- sd(daily$amp, na.rm=T)
  
}


for(i in 1:9){
  site <- NHCsites2018[i]
  daytimes <- data.frame(DO=DOdat[[paste(site,"DO_mgL", sep=".")]],
                      hour=DOdat$hour, datetime = DOdat$DateTime_UTC)
  daytimes <- daytimes[-which(daytimes$hour>16),]
  daytimes <- daytimes[-which(daytimes$hour<9),]
  
  DOSiteSummary$dayMean[i] <- mean(daytimes$DO, na.rm=T)
}

#save summary file
write.csv(DOSiteSummary, file = "data/DOtimeseriesMetrics.csv", row.names=FALSE)


# Model ER using nighttime regression from Stream Metabolizer
model_name <- mm_name(type='night')
model_specs <- specs(model_name)
mm <- metab(model_specs, data=DOdat[,c(1,2)])


write.csv(DOSiteSummary, file="allsites_DODailySum.csv")
par(mfrow = c(1,1))
par(mar = c(3,4,1,3))
barplot(DOSiteSummary$PrHypox5, col = "grey80", ylab = "Pr[DO < x]",
        names.arg = NHCsites2018)#, ylim = c(0,0.85))
par(new = T)
barplot(DOSiteSummary$PrHypox2, col = "grey40", axes=FALSE,
        ylim = c(0,0.85), bty = 'n')
par(new = T)
barplot(DOSiteSummary$PrAnox, col = "grey10", axes=FALSE, 
        ylim = c(0,0.85),bty = 'n')
legend("topright", bty = 'n', legend = c("Pr[DO<5]","Pr[DO<2]","Pr[DO<0]"),
       col  = c("grey80", "grey40","grey10"), pch = 15) 

plot(DOSiteSummary$cross0, ylab = "Return Interval (days)", ylim = c(0,45),
         type = "b",names.arg = NHCsites2018, pch = 1, lwd = 2,col = "black")
lines(DOSiteSummary$cross2, col = "grey40", lwd = 2, type = "b")
lines(DOSiteSummary$cross5, col = "grey80", lwd = 2, type = "b")
legend("topleft", bty = 'n', legend = c("5 mgL","2 mgL","0 mgL"),
       col  = c("grey80", "grey40","black"), pch = 1, lwd = 3) 



# Data from MC751
DOSiteSummary[which(DOSiteSummary>365)]
site <- "MC751"
site_code <- paste('NC',site,sep="_")

# Get list of variables and start and end date for each site from the database
mdat <- query_available_data(region = 'NC',site = site)
variables <- as.vector(mdat$variables$variables)
start_date <- mdat$datebounds$firstRecord
end_date <- mdat$datebounds$lastRecord

dat <- request_data(sitecode=site_code, variables = variables,
                    startdate=start_date, enddate=end_date)
w <- which(dat$data$value< -0.5)
dat$data$value[w]<-NA
dat <- dat$data[,c(1,4,5)]%>% spread(key = "variable",
                                          value = "value")
dat.xts <- xts(dat[,c(2,3,8,9)], order.by = dat[,1])
dygraph(dat.xts[,1])
