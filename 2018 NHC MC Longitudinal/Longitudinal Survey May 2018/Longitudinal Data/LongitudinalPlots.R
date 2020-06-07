### Code to read in and plot data from longitudinal DO sampling in May 2018
# Created by Alice Carter
# Nov 12 2018

# Significant Updates: Oct 28 2018
#          - finalized dataset as: LongitudinalSamples_2018May.csv


library(tidyr)
library(lubridate)
#library(RColorBrewer)
#library(SDMTools)
library(dplyr)
library(zoo)

# Load Longitudinal summary datafile
# Called CompiledLongitudinalSamples_May2018.csv, it is a sheet in 2018May_AliceNetworkSampling.xlsx
dat <- read.csv(file = "LongitudinalSamples_2018May.csv", header = T, 
                stringsAsFactors = FALSE, na.strings = c("", "NA"))
dat<- dplyr::rename(dat, streamSection=colnames(dat)[1])

# Make sure dates and times are properly formatted
dat$datetime<- mdy_hm(paste(dat$Date,dat$Time, sep=" "))
dat$Time<- as.POSIXct(dat$Time, format = "%H:%M")
dat$Date <- mdy(dat$Date)


colors <- c("#2264A8", "#6BA5D7","#CBDBE5","#8D9295","#F56146","#FFAC88","#1B325F")


#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c(colors[7],"white"))

#This adds a column of color values
# based on the y collection times
dat$color <- rbPal(10)[as.numeric(cut(dat$Time,breaks = 10))]
minTime <- as.POSIXct("07:23", format="%H:%M")
maxTime <- as.POSIXct("19:00", format="%H:%M")

# Subset dat based on which sections you want to plot.
# MC_trib follows the path up to the headwaters of mud tributary while MC_751 and MCU_ConfAmVi follow a trib through american village

dat751 <- filter(dat, streamSection != "MC_trib")
datMtrib <- filter(dat, streamSection != "MCU_ConfAmVi") %>% filter(streamSection != "MC_751")

dat <- datMtrib

roads <- filter(dat, Road.Crossing == 1) %>% select(distance_m)
wwtp <- filter(dat, WWTP ==1) %>% select(distance_m)
snsrs <- filter(dat, !is.na(SampleStation)) %>% select( SampleStation, distance_m, DO_mgL, DO_pctsat)

# prep nutrient data
dat.ions <- dat %>% select(datetime, distance_m, SampleStation, WWTP, DO_pctsat, DO_mgL, SpC_uScm, Cl.mgL, NO3.N.mgL, NH4.N.mgL, SO4.mgL)
dat.ions <- dat.ions[order(dat.ions$distance_m),]
dat.ions$distance_m <- dat.ions$distance_m/1000
dat.upper<- dat.ions[which(!is.na(dat.ions$Cl.mgL)),]
dat.upper<- rbind(dat.ions[1,], dat.upper)
w <- which(dat.upper$distance_m >= wwtp$distance_m/1000)
dat.lower<- dat.upper[w,]
dat.upper[w,8:11]<- NA

dat.sum <- colMeans(dat.lower[,c(2,9,10,11)], na.rm=T)


png(file="LongDOChemplot2.png", width=8, height=4.8, units="in", type="cairo", res=300)
  
  par(mar =c(0,0,0,5), mgp=c(2,1,0), oma=c(6,4,1,0), mfrow = c(2,1))
  
  
  plot(dat$distance_m/1000, dat$DO_pctsat, pch = 21,cex=0.7,yaxt="n",xaxt="n", ylim = c(0,140))
  axis(2, at=c(0,25,50,75,100), cex.lab=.9,cex.axis =.9)
  polygon(x = c(-2,30,30,-2),y = c(-10,-10,50,50), col = "grey90", border = NA )
  #polygon(x = c(-2,30,30,-2),y = c(-2,-2,3,3), col = "grey70", border = NA )
  #abline(v = roads$distance_m/1000, col =' black',lwd = 1.5, cex = 4, lty = 2)
  abline(v=wwtp$distance_m/1000, col = "steelblue", lwd = 2, lty = 1)
  par(new="T")
  plot(dat$distance_m/1000, dat$DO_pctsat, ylim = c(0,140),cex=0.7,pch = 21, bg = dat$color,xaxt='n', yaxt='n', ann=F)
  points(snsrs$distance_m/1000, snsrs$DO_pctsat, pch = 20, cex = 1.5, col = colors[5])
  #text(x=wwtp$distance_m/1000, y=3, labels="WWTP", pos=4, offset=0.4, cex=.8)
  
  gradientLegend(valRange=c(0,1),
                 color = rbPal(10), pos = c(21.5, 5,22.2,135), side = 4,length = 0.5, depth = 0.03, inside = FALSE, 
                 coords = TRUE, pos.num = NULL, n.seg =0, border.col = "black", dec = NULL,
                 fit.margin = FALSE)
  text(x=22.2, y=7, labels="7:00", pos=4, offset=0.2, xpd=NA, cex=.7)
  text(x=22.2, y=131, labels="19:00", pos=4, offset=0.1, xpd=NA, cex=.7)
  mtext(text = "sample time", side=4, line=1.5, outer=FALSE, cex=.7)
  mtext("DO (%sat)", 2, line=2)
  legend(0,140, legend=c("DO (%sat)    ","sensor  ","WWTP", "hypoxic"), 
         cex=.8, bty="n",ncol=4, xpd=NA, 
         pt.cex=c(.8,1),
         border=NA,
         col=c("black",colors[5],"steelblue",NA),
         lty=c(NA,NA,1,NA),
         pch=c(1,19,NA,NA),
         fill=c(NA,NA,NA,"grey90"))
  
  
  ################################################################################
  # Nutrient Concentrations normalized to Cl and relative to mean
  
  DO.col <-rbPal(3)[2]
  NH4.col <- 3
  NO3.col <- 2
  SO4.col <-1
  
  plot(dat.upper$distance_m, dat.upper$SO4.mgL, type="l", lty=SO4.col, lwd=2,
       xlim= c(0,max(dat.upper$distance_m)), ylim = c(0,10), xlab="", ylab="", xaxt="n", yaxt="n")
  polygon(c(dat.upper$distance_m, rev(dat.upper$distance_m)), na.approx(c(dat.upper$DO_pctsat/20, rep(0, nrow(dat.upper)))),
          border=F, col = alpha(DO.col,.4))
  lines(dat.upper$distance_m, dat.upper$NO3.N.mgL, lwd=2,lty=NO3.col)
  axis(2, at=c(0,2,4,6,8), cex.lab=.9,cex.axis =.9)
  mtext("NO3-N, SO4, (mg/L)", side=2, line=2)
  
  text(19.2, 4.5, "NH4 = 0.05 mg/L",  cex=.8, pos=1)
  text(19.2, 5.5, "NO3 = 29 mg/L", cex=.8, pos=1)
  text(19.2, 6.5, "SO4 = 61 mg/L", cex=.8, pos=1)
  text(19.2, 9, "mean concentration\nbelow WWTP", cex=.7)
  
  par(new=T)
  
  plot(dat.upper$distance_m, dat.upper$NH4.N.mgL, type="l", lwd=2, lty =NH4.col,
       xaxt="n", yaxt="n", ylab="", xlab="", ylim = c(0,.2))
  axis(4, at=c(0,.06,.12,.18), cex.lab=.9,cex.axis =.9)
  
  abline(v=wwtp$distance_m/1000, col="steelblue", lwd=2)
  mtext("NH4-N, (mg/L)", side=4, line=2)
  axis(1, at=c(0,4,8,12,16,20), cex.lab=.9)
  mtext("Distance downstream (km)", side=1, line=2)
  legend(0,.2, legend=c("DO (%sat)","SO4","NO3-N", "NH4-N"), 
         cex=.8, bty="n",ncol=4, xpd=NA, 
         pt.cex=1.2,
         fill =c(DO.col,NA,NA,NA),
         border=NA,
         col=c("white","black","black","black"),
         lty=c(NA,1,2,3))
  

    
  dev.off()
##################################################################################
# Plot Nutrient Data
# Full Nutrient Concentration data

par(mar = c(4,4,1,1), oma = c(1,1,0,0))

cols <- c("orange", "purple", "darkred", "black")
plot(dat$distance_m/1000, dat$SO4.mgL, pch = 20,
      ylab = "concentration (mg/L)", xlab = "Kilometers downstream", log="y")
# w2 <- which(dat$DO_mgL<2)
# w5 <- which(dat$DO_mgL<5)

# abline(v=dat$distance_m[w5]/1000, col = "grey90", lwd = 3)
# abline(v=dat$distance_m[w2]/1000, col = "grey60", lwd = 3)
par(new=T)

plot(dat$distance_m/1000, dat$SO4.mgL, pch = 20, col = cols[2], cex=1.2,
     xaxt="n", yaxt="n",ann=F, ylim = c(.1, 55),log="y")

points(dat$distance_m/1000, dat$NO3.N.mgL, pch = 20, col = cols[1], cex=1.2)
points(dat$distance_m/1000, dat$NH4.N.mgL, col =cols[3], pch = 20, cex=1.2)
points(dat$distance_m/1000, dat$Cl.mgL, col = cols[4], pch = 20, cex=1.2)
legend("topleft", c("NO3", "SO4","NH4", "Cl"), pch=c(20,20,20,20), col = cols, bg='n', bty="n")



dat.norm <- dat.upper%>% select(distance_m, 
                          NO3Clnorm=NO3.N.mgL, 
                          NH4norm=NH4.N.mgL, Cl.mgL, 
                          SO4norm=SO4.mgL, DO_pctsat)
dat.norm$NO3Clnorm <- dat.norm$NO3Clnorm/dat.norm$Cl.mgL
dat.norm$NO3Clnorm <- dat.norm$NO3Clnorm/mean(dat.norm$NO3Clnorm, na.rm=T)

dat.norm$NH4norm <- dat.norm$NH4norm/dat.norm$Cl.mgL
dat.norm$NH4norm <- dat.norm$NH4norm/mean(dat.norm$NH4norm, na.rm=T)
dat.norm$SO4norm <- dat.norm$SO4norm/dat.norm$Cl.mgL
dat.norm$SO4norm <- dat.norm$SO4norm/mean(dat.norm$SO4norm, na.rm=T)


dat.norm$Cl.mgL <- dat.norm$Cl.mgL/mean(dat.norm$Cl.mgL, na.rm=T)
dat.norm <- dat.norm[which(!is.na(dat.norm$Cl.mgL)),]
dat.norm <- dat.norm[order(dat.norm$distance_m),]
ylims <- c(0.1,5)
plot(dat.norm$distance_m, dat.norm$NO3Clnorm, type="l",ylim=c(0,2), xlim = c(1000,16200))#, log="y")
lines(dat.norm$distance_m, dat.norm$SO4norm, col=2)
lines(dat.norm$distance_m, dat.norm$NO3Clnorm, col=3)
lines(dat.norm$distance_m, dat.norm$NH4norm, col=4)
polygon(c(dat.norm$distance_m, rev(dat.norm$distance_m)), 
        na.approx(c(dat.norm$DO_pctsat/50, rep(0,nrow(dat.norm)))),
        col = alpha(1, .3), border=F)
#################################################################################
# relative to mean, excluding points below the wwtp
# Remove data below the wwtp
dat.ions <- dat %>% select(datetime, distance_m, SampleStation, WWTP, DO_pctsat, DO_mgL, SpC_uScm, Cl.mgL, NO3.N.mgL, NH4.N.mgL, SO4.mgL)
w <- which(dat$distance_m >= wwtp$distance_m)
dat.ions[w,8:11] <- "NA"
m<- max(as.numeric(dat.ions[8:11]), na.rm=T)
par(mar = c(4,4,1,5), oma = c(1,1,0,0))

cols <- c("orange", "purple", "darkred", "black")
plot(dat.ions$distance_m/1000, dat.ions$SO4.mgL, pch = 20,ylim = c(0,8),
     ylab = "NO3, SO4 (mg/L)", xlab = "Kilometers downstream")
w2 <- which(dat.ions$DO_mgL<2)
w5 <- which(dat.ions$DO_mgL<5)

abline(v=dat.ions$distance_m[w5]/1000, col = "grey90", lwd = 3)
abline(v=dat.ions$distance_m[w2]/1000, col = "grey60", lwd = 3)
par(new=T)

plot(dat.ions$distance_m/1000, dat.ions$SO4.mgL, pch = 20, col = cols[2], cex=1.2,
     xaxt="n", yaxt="n",ann=F, ylim = c(0,8))

points(dat.ions$distance_m/1000, dat.ions$NO3.N.mgL, pch = 20, col = cols[1], cex=1.2)

par(new = T)
plot(dat.ions$distance_m/1000, dat.ions$NH4.N.mgL, col =cols[3], pch = 20, cex=1.2, xaxt='n', yaxt="n", ann=F)
axis(4)
mtext("NH4 (mg/L)", side=4, line=3, cex.lab=1, las=3)

legend("topright", c("NO3", "SO4","NH4"), pch=c(20,20,20), col = cols[1:3], bg='n', bty="n")




dat.ions$Cl.dev <- as.numeric(dat.ions$Cl.mgL)/mean(as.numeric(dat.ions$Cl.mgL), na.rm=T)
dat.ions$SO4.dev <- as.numeric(dat.ions$SO4.mgL)/mean(as.numeric(dat.ions$SO4.mgL), na.rm=T)
dat.ions$NH4.dev <- as.numeric(dat.ions$NH4.N.mgL)/mean(as.numeric(dat.ions$NH4.N.mgL), na.rm=T)
dat.ions$NO3.dev <- as.numeric(dat.ions$NO3.N.mgL)/mean(as.numeric(dat.ions$NO3.N.mgL), na.rm=T)

par(mar = c(4,4,1,1), oma = c(1,1,0,0))

cols <- c("orange", "purple", "darkred", "black")
m <- max(dat.ions[12:15], na.rm=T)
plot(dat.ions$distance_m/1000, dat.ions$SO4.dev, pch = 20,ylim = c(0,m),
     ylab = "concentration (mg/L)", xlab = "Kilometers downstream")
w2 <- which(dat$DO_mgL<2)
w5 <- which(dat$DO_mgL<5)

abline(v=dat$distance_m[w5]/1000, col = "grey90", lwd = 3)
abline(v=dat$distance_m[w2]/1000, col = "grey60", lwd = 3)
par(new=T)

plot(dat.ions$distance_m/1000, dat.ions$SO4.dev, pch = 20, col = cols[2], cex=1.2,
     ylim = c(0,m),xaxt="n", yaxt="n",ann=F)

points(dat.ions$distance_m/1000, dat.ions$NO3.dev, pch = 20, col = cols[1], cex=1.2)
points(dat.ions$distance_m/1000, dat.ions$NH4.dev, col =cols[3], pch = 20, cex=1.2)
points(dat.ions$distance_m/1000, dat.ions$Cl.dev, col = cols[4], pch = 20, cex=1.2)









plot(dat$distance_m/1000, dat$NH4_mgL/dat$NO3_mgL, ylim = c(0,0.2), pch = 20,
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
