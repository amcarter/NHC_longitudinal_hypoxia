### Code to read in and plot data from longitudinal DO sampling in May 2018
# Created by Alice Carter
# Nov 12 2018

# Significant Updates: Oct 28 2018
#          - finalized dataset as: LongitudinalSamples_2018May.csv


library(tidyr)
library(RColorBrewer)
library(SDMTools)
library(dplyr)

# Load Longitudinal summary datafile
# Called CompiledLongitudinalSamples_May2018.csv, it is a sheet in 2018May_AliceNetworkSampling.xlsx
dat <- read.csv(file = "LongitudinalSamples_2018May.csv", header = T, 
                stringsAsFactors = FALSE, na.strings = c("", "NA"))


#colors <- c(,"#9CC4E4","#E9F2F9","#3A89C9","#F26C4F")
colors <- c("#2264A8", "#6BA5D7","#CBDBE5","#8D9295","#F56146","#FFAC88","#1B325F")

roads <- filter(dat, Road.Crossing == 1) %>% select(distance_m)
wwtp <- filter(dat, WWTP ==1) %>% select(distance_m)
snsrs <- filter(dat, !is.na(SampleStation)) %>% select( SampleStation, distance_m)

#cols <- brewer.pal(3, "Dark2")
#shades <- c("#bbc7e0",cols[3],"#5e79b6")

#hypox5 <- nrow(filter(dat,DO_mgL<5))/nrow(dat)
#hypox2 <- nrow(filter(dat,DO_mgL<2))/nrow(dat)


#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c(colors[7],"white"))

#This adds a column of color values
# based on the y collection times
dat$Col <- rbPal(10)[as.numeric(cut(dat$Time,breaks = 10))]





par(mar =c(4,5,2,1))
plot(dat$distance_m/1000, dat$DO_mgL, pch = 20, cex.lab = 1.8,cex.axis = 1.15,
     ylab = "DO (mg/L)", xlab = "Kilometers downstream")
polygon(x = c(-2,30,30,-2),y = c(-2,-2,6,6), col = "grey90", border = NA )
polygon(x = c(-2,30,30,-2),y = c(-2,-2,3,3), col = colors[1], border = NA )

points(dat$distance_m/1000, dat$DO_mgL, pch = 20)
abline(v = roads$distance_m/1000, col = colors[1],lwd = 1.7, cex = 4)
abline(v=wwtp$distance_m/1000, col = colors[5], lwd = 2, lty = 2)
points(snsrs$distance_m/1000, snsrs$SampleStation, pch = 20, cex = 2.5, col = colors[5])
#text(x = snsrs$distance_m/1000, y = c(9,9,9), labels = c("A","B","C"))
legend("topleft", c("major road crossing", "waste water"), cex = 1.2,lty = c(1,2), lwd = 2, col = colors[c(1,5)])

legend("topleft", c("sensor location"), pch = c(19),col = colors[5])
legend("topleft", c("7:00","19:00", "Time"), pch = c(19,19,19),col = colors[5])

# Plot Nutrient Data
par(mar = c(4,4,1,1), oma = c(1,1,0,0))
plot(dat$distance_m/1000, (dat$DO_mgL), pch = 20, cex.lab = 1,cex.axis = 1,
     col = "grey80")#, ylim = c(0,32))
plot(dat$distance_m/1000, (dat$NO3_mgL/dat$Cl_mgL)*30, pch = 20, col = "orange",
     ylim = c(0,10),  ylab = "concentration", xlab = "Kilometers downstream")
w2 <- which(dat$DO_mgL<2)
w5 <- which(dat$DO_mgL<5)

abline(v=dat$distance_m[w5]/1000, col = "grey80", lwd = 2)
abline(v=dat$distance_m[w2]/1000, col = "grey50", lwd = 2)

points(dat$distance_m/1000, (dat$NO3_mgL/dat$Cl_mgL)*30, pch = 20, col = "orange")
points(dat$distance_m/1000, (dat$SO4_mgL), pch = 20, col = "purple", ylim = c(0,60))
points(dat$distance_m/1000, dat$NH4_mgL*30, col = "darkred", pch = 20)
legend("topright", c("NO3", "SO4","NH4"), pch=c(20,20,20), col = c("orange","purple","darkred"), bty="n")


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
