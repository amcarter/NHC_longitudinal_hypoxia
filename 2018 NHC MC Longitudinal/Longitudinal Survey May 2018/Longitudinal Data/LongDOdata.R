### Code to read in and plot data from longitudinal DO sampling in May 2018
# Created by Alice Carter
# Nov 12 2018

# Significant Updates: Oct 28 2018
#          - finalized dataset as: LongitudinalSamples_2018May.csv

library(tidyr)
library(dplyr)
library(readr)
library(lubridate)
library(streamMetabolizer)
setwd(hypox_projdir)

# Load Longitudinal summary datafile
# Called CompiledLongitudinalSamples_May2018.csv, it is a sheet in 2018May_AliceNetworkSampling.xlsx
dat <- read_csv(file = "data/long_survey_withNHD.csv")


roads <- filter(dat, RoadCrossing == 1) %>% select(distance_m)
wwtp <- filter(dat, WWTP ==1) %>% select(distance_m)
snsrs <- filter(dat, !is.na(SampleStation)) %>% select( SampleStation, distance_m)

#cols <- brewer.pal(3, "Dark2")
#shades <- c("#bbc7e0",cols[3],"#5e79b6")

covars <- dat %>%
  select(Time, streamSection,distance_m, slope, Latitude, Longitude, width_m, depth_m, velocity_ms, temp_C,
         DO_pctsat, DO_mgL, Cl.mgL, SO4.mgL, NO3.N.mgL, NH4.N.mgL, Habitat)

sunrise <- hms("06:04:00")

covars$light_time <- hms(as.character(covars$Time))-sunrise
covars$light_hrs <- hour(covars$light_time)+minute(covars$light_time)/60
covars$Habitat[1:12]<- "Ri"
covars$Habitat <- as.factor(covars$Habitat)
boxplot(DO_pctsat~Habitat, data=covars, plot=F)

#######################################
#rescale predictor variables
covars$slope_mkm <- covars$slope*1000
covars$DO.upstream <- c(covars$DO_pctsat[1], covars$DO_pctsat[1:(nrow(covars)-1)])

mm <- lmer(DO_pctsat~DO.upstream+velocity_ms+ slope_mkm +(1|streamSection), data=covars)
confint(mm)
summary(mm)
AIC(mm)
covars$DOpred <- 11.3085+.6971*covars$DO.upstream+42.6263*covars$velocity_ms+1.7101*covars$slope_mkm

mm1 <- lmer(DO_pctsat~light_hrs+DO.upstream+velocity_ms+ slope_mkm +(1|streamSection), data=covars)
confint(mm1)
summary(mm1)
AIC(mm1)
covars$DOpred <- 7.66648+.62363*covars$DO.upstream+40.27201*covars$velocity_ms+1.45674*covars$slope_mkm+1.36948*covars$light_hrs
covars$pred.lower <- -.512383+.5199866*covars$DO.upstream+18.3703455*covars$velocity_ms+.4244931*covars$slope_mkm+.5155646*covars$light_hrs
covars$pred.upper <- 16.1084988+.7511535*covars$DO.upstream+62.55656*covars$velocity_ms+2.3323187*covars$slope_mkm+2.1859767*covars$light_hrs
covars <- covars[order(covars$distance_m),]

png("figures/long_DO_model.png", width=8, height = 5, units="in", res=300)
par(mar=c(4,4,1,1), oma=c(0,0,2,0)) 
plot(covars$distance_m/1000, covars$pred.upper, ylim = c(0,128),type="n" , xlab="distance (km)", ylab = "DO (% sat)")
lines(covars$distance_m/1000, covars$DOpred, lwd=2, col="steelblue")
polygon(c(covars$distance_m/1000, rev(covars$distance_m/1000)), 
        na.approx(c(covars$pred.lower, rev(covars$pred.upper)),na.rm=F), 
        col=alpha("steelblue",.3), border=NA)
points(covars$distance_m/1000, covars$DO_pctsat, pch=20)
par(new=T, oma=c(0,0,0,0))
legend("top",legend=c("Measured DO","Modeled DO", "95% CI"),
       fill=c(NA,NA, alpha("steelblue",.3)),
       col=c(1,"steelblue",NA),xpd=NA,
       border=NA,bty="n",ncol=3,
       lty=c(NA,1,NA),
       lwd=c(NA,2,NA),
       pch=c(20,NA,NA))

dev.off()
covars$NH4.N.ugL<- covars$NH4.N.mgL*1000
no3 <- lm(NO3.N.mgL~DO_pctsat, data=covars)
nh4 <- lm(NH4.N.ugL~DO_pctsat, data=covars)



plot(covars$DO_pctsat- covars$DOpred)
points(covars$DO_pctsat- covars$DOpred1, col=2)

plot(dat$DO_pctsat, dat$NH4.N.mgL/dat$NO3.N.mgL, ylim = c(0,0.2), pch = 20,
     ylab = "NH4/NO3")

abline(v=dat$distance_m[w5]/1000, col = "grey80", lwd = 2)
abline(v=dat$distance_m[w2]/1000, col = "grey50", lwd = 2)
points(dat$distance_m/1000, dat$NH4_mgL/dat$NO3_mgL, ylim = c(0,0.2), pch = 20)

plot(dat$distance_m/1000, dat$NH4_mgL/dat$NO3_mgL, ylim = c(0,0.2), pch = 20, ylab = "CH4/CO2")
abline(v=dat$distance_m[w5]/1000, col = "grey80", lwd = 2)
abline(v=dat$distance_m[w2]/1000, col = "grey50", lwd = 2)
points(dat$distance_m/1000, dat$NH4_mgL/dat$NO3_mgL, ylim = c(0,0.2), pch = 20)


axis(side=4)
mtext(side = 4, line = 3, 'NO3, SO4 (mg/L)')
par(new=T)
points(dat$distance_m/1000, (dat$NH4_mgL), pch = 20,xlab = NA, ylab = NA, axes=F, col = "darkred")
points(dat$distance_m/1000, dat$Cl_mgL/2, pch = 20, col = "black")
par(new = T)
plot(dat$distance_m/1000, dat$Cl_mgL, ylim = c(10,20))
mtext(side = 4, line = 3, 'NH4 (mg/L)')
legend("topleft", col = c("grey80", "orange", "darkred","purple", "black"), pch = 20, legend = c("O2","NO3","NH4","SO4", "Cl"))

plot(dat$DO_mgL, log(dat$NO3_mgL/dat$NH4_mgL))
plot(dat$DO_mgL, (dat$SO4_mgL/dat$Cl_mgL))

plot(dat$DO_mgL, dat$SO4_mgL/dat$Cl_mgL)#, ylim = c(0.1,0.5))
plot(dat$SO4_mgL/dat$Cl_mgL, dat$NO3_mgL/dat$NH4_mgL, xlab = "SO4/Cl",
     ylab = "NO3/NH4", pch = 19
     ,ylim = c(0,120), xlim = c(0.1,0.5) )

plot(dat$DO_mgL, dat$SO4_mgL, col = "purple", pch = 19, xlab = "DO mgL",
     ylab = "SO4 mgL", ylim = c(1,9))






# Water Quality Portal Dat

EPA <- read_csv("DOpiedmontEPAWQ.csv")

EPA$Date <- as.POSIXlt(EPA$Date, format = "%m/%d/%y")

EPA$Date[which(EPA$Date >as.POSIXlt("2060-01-01"))]$year <-
  EPA$Date[which(EPA$Date >as.POSIXlt("2060-01-01"))]$year - 100

EPA <- EPA[-which(EPA$DO>(mean(EPA$DO)+4*sd(EPA$DO))),]
plot(EPA$Date, EPA$DO)


EPA$hypox5 <- rep(0, nrow(EPA))
EPA$hypox5[which(EPA$DO<5)] <- 1

EPA$hypox2 <- rep(0, nrow(EPA))
EPA$hypox2[which(EPA$DO<2)] <- 1

EPA$month <- EPA$Date$mon
EPA$Year <- EPA$Date$year
EPA$Date <- as.Date(EPA$Date)
Jun <- filter(EPA, month==6)
Aug <- filter(EPA, month == 8)
Jul <- filter(EPA, month==7)
Sept <- filter(EPA, month == 9)
Oct <- filter(EPA, month == 10)

mons <- c("Jul","Aug","Sept","Oct")
hypox <- data.frame(month = mons, hypox2 = NA, hypox5 = NA)
tmp <- Jun
hypox[3,2] <- sum(tmp$hypox2)/nrow(tmp)
hypox[3,3] <- sum(tmp$hypox5)/nrow(tmp)


HypoxMonth <- EPA %>%
  group_by(month) %>%
  summarise(hypox2 = mean(hypox2, na.rm = T),
            hypox5 = mean(hypox5, na.rm = T),
                   hypox2.stdev = sd(hypox2, na.rm = T),
                   hypox5.stdev = sd(hypox5, na.rm = T))

HypoxYear <- EPA %>%
  group_by(Year)%>%
  dplyr::summarize(hypox2 = mean(hypox2, na.rm = T),
                   hypox5 = mean(hypox5, na.rm = T),
                   hypox2.stdev = sd(hypox2),
                   hypox5.stdev = sd(hypox5))
HypoxYear <- EPA %>%
    group_by(Year) %>%
    summarise_each(funs(mean,sd))
HypoxYear$Year <- HypoxYear$Year +1900
HypoxMonth <- EPA %>%
    group_by(month) %>%
    summarise_each(funs(mean,sd))
plot(HypoxMonth$month, HypoxMonth$hypox5_mean, xlab = "Month", ylab = "Pr[hypoxia]",
     ylim = c(0,max(HypoxMonth$hypox5_mean+HypoxMonth$hypox5_sd)))
polygon(HypoxMonth$month, (HypoxMonth$hypox5_mean - HypoxMonth$hypox5_sd),
       HypoxMonth$month, (HypoxMonth$hypox5_mean + HypoxMonth$hypox5_sd),
        col = "lightgrey")
points(HypoxMonth$month, HypoxMonth$hypox2_mean, pch = 19)
arrows(HypoxMonth$month, (HypoxMonth$hypox2_mean - HypoxMonth$hypox2_sd),
       HypoxMonth$month, (HypoxMonth$hypox2_mean + HypoxMonth$hypox2_sd),
       length = 0.1, angle = 90)

legend("topright", legend = c("below 5 mg/L","below 2 mg/L"), pch = c(1,19))



plot(HypoxMonth$month, HypoxMonth$hypox5_mean, xlab = "month", ylab = "Pr[hypoxia]")
points(HypoxMonth$month, HypoxMonth$hypox2_mean, pch = 19)
legend("topright", legend = c("below 5 mg/L","below 2 mg/L"), pch = c(1,19))



plot(HypoxYear$Year, HypoxYear$hypox5_mean, xlab = "Date", ylab = "Pr[hypoxia]",
     xlim = c(1968,2016))
points(HypoxYear$Year, HypoxYear$hypox2_mean, pch = 19)
legend("topright", legend = c("below 5 mg/L","below 2 mg/L"), pch = c(1,19))


EPA$Year<-  as.nueric(substr(EPA$Date, 1,4))
plot(EPA$Date, EPA$Time)#, col = EPA$Time)e



#############
# Read NHC data and plot

dat <- read.csv("NC_NHC_sensorData.csv")
LAI <- read.csv("NC_UNHC_2017.csv")
dat[which(dat$flagtype=="Bad Data"),5] <- NA
dat$DateTime_UTC <- as.POSIXct(dat$DateTime_UTC, format = "%Y-%m-%d %H:%M:%S")
dat <- dat[-which(is.na(dat$DateTime_UTC)),]
dat <- spread(dat, variable, value)
dat <- dat[,c(-2,-3,-4,-5)]

dat$lvl_m <- 0.101972*(dat$WaterPres_kPa-dat$AirPres_kPa)+0.24
dat$lvl_m[which(dat$lvl_m <0.3)]<- NA

dat$disc_cms <- 0.4386*(dat$lvl_m^7.7976)

dat$DO_mgL[which(dat$DO_mgL>20)] <- NA
daterng <- c(min(dat$DateTime_UTC, na.rm = T)-50000000,max(dat$DateTime_UTC, na.rm = T)+10000000)

par(mfrow = c(1,1), mar = c(4,5,1,3))

plot(dat$DateTime_UTC, dat$DO_mgL, pch = 20, cex = .7, ylab = "DO (mg/L)", cex.axis = 1.3,
     cex.lab = 1.3, ylim = c(0,15), xlab = "LAI")
polygon(x = c(daterng[1],daterng[2],daterng[2],daterng[1]),y = c(-2,-2,5,5), col = "grey90", border = NA )
points(dat$DateTime_UTC, dat$DO_mgL, pch = 20, cex = .7)
par(new = TRUE)
plot(LAI$LAI_proc, type = "l", axes = FALSE, col = "darkgreen", lwd = 2.5, ylim = c(0,5.5),ylab = "")
axis(4, at = seq(0,5, by = 1), cex.axis = 1.3)
plot(dat$DateTime_UTC, dat$disc_cms, type = "l", lwd = 2, col = colors[7],
     xlim = dates, log = "y",ylab = "Discharge\n(cms)", ylim = c(.001,4000))
plot(dat$DateTime_UTC, dat$WaterPres_kPa, pch = 20)

plot(dat$DateTime_UTC, dat$AirPres_kPa, type = "l", lwd = 2, col = colors[7],
     xlim = dates,ylab = "Discharge\n(cms)")

par(mar = c(2,4,1,1),mfrow = c(3,1))

dates <- as.POSIXct(c("2017-09-20 00:00:00","2017-10-20 00:00:00"))
plot(dat$DateTime_UTC, dat$DO_mgL, pch = 20, xlim = dates,cex = .8,
     cex.lab = 1.3, cex.axis = 1.3,ylab = "DO (mg/L)", ylim = c(0,15))
axis(side = 1, at = seq(dates[1],to = dates[2],length.out = 4), labels = c("Sep-20","Sep-30","Oct-10","Oct-20"))
polygon(x = c(daterng[1],daterng[2],daterng[2],daterng[1]),y = c(-2,-2,5,5), col = "grey90", border = NA )
points(dat$DateTime_UTC, dat$DO_mgL, pch = 20, cex = .7)





boxplot(HypoxYear[,2:3]*100,ylab = "Percent of Observations per Year", col = c("grey30","grey90"))
legend("topleft", legend = c("Less than 2 mg/L","Less than 5 mg/L"))
