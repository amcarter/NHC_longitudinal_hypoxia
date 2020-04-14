### Code to generate semi-variogram from longitudinal DO sampling in May 2018
# Created by Alice Carter
# Jan 15 2020

library(tidyverse)
library(lubridate)
library(gstat)

# Load Longitudinal summary datafile
# Called CompiledLongitudinalSamples_May2018.csv, it is a sheet in 2018May_AliceNetworkSampling.xlsx
dat <- read.csv(file = "LongitudinalSamples_2018May.csv", header = T, 
                stringsAsFactors = FALSE, na.strings = c("", "NA"))

# Make sure dates and times are properly formatted
dat$datetime<- mdy_hm(paste(dat$Date,dat$Time, sep=" "))
glimpse(dat)
# build distance matrix
dat$x <- dat$distance_m/1000
dat$y <- 0


# Spatialdat <- dat %>% select(x,y, DO_mgL)
# w <- which(is.na(Spatialdat$DO_mgL))
# Spatialdat <- Spatialdat[-w,]
# wwtp <- filter(dat, WWTP ==1)$distance_m
# 
# aboveWWTP<- Spatialdat[which(Spatialdat$x<wwtp/1000),]

# Create vector of NHC and MC data separately
unique(dat$ï..streamSection)
MC <- filter(dat, ï..streamSection %in% c("MC_trib","MC_751","MC_Erwin","MC_Pickett"))%>%
  select(x,y,DO_mgL) %>% mutate(x=x-min(x, na.rm=T))
MC <- MC[which(!is.na(MC$DO_mgL)),]
NHC <- filter(dat, ï..streamSection %in% c("NHC_University", "NHC_KingCharles","NHC_Stagecoach"))%>%
  select(x,y,DO_mgL)%>% mutate(x=x-min(x, na.rm=T))
NHC <- NHC[which(!is.na(NHC$DO_mgL)),]






# basic variogram
MC.vgm <- variogram(DO_mgL~1, loc= ~x+y, data=MC, width = 0.1)
MC.fit <- fit.variogram(MC.vgm, vgm("Gau"))
GauModel <- vgm(psill=MC.fit$psill[2], model="Gau", nugget=MC.fit$psill[1], range=MC.fit$range[2])
ExpModel <- vgm(psill=7.0136447, model = "Exp", nugget = 0.4112233, range = 1.977862)
plot(MC.vgm, model = GauModel, col = "black", pch = 19, ylab = "gamma", xlab = "distance (km)")

NHC.vgm <- variogram(DO_mgL~1, loc= ~x+y, data=NHC, width = 0.1)
NHC.fit <- fit.variogram(NHC.vgm, vgm("Gau"))
GauModel <- vgm(psill=NHC.fit$psill[2], model="Gau", nugget=NHC.fit$psill[1], range=NHC.fit$range[2])
plot(NHC.vgm, model = GauModel, col = "black", pch = 19, ylab = "gamma", xlab = "distance (km)", ylim = c(0,4))



DO.vgm1 <- variogram(DO_mgL~1, loc= ~x+y, data=aboveWWTP)
DO.fit1 <- fit.variogram(DO.vgm1, vgm("Gau","Sph","Exp"))
GauModel <- vgm(psill=5.6182653, model="Gau", nugget=0.7965, range=1.2622)
plot(DO.vgm1, model = GauModel, col = "black", pch = 19, ylab = "gamma", xlab = "distance (km)")

# Conductivity variogram
SpcDat <- dat[-which(is.na(dat$SpC_uScm)),]
SpcDat <- SpcDat[which(SpcDat$x<wwtp/1000),]
Spc.vgm <- variogram(SpC_uScm~1, loc = ~x+y, data = SpcDat, width = .1)
Spc.fit <- fit.variogram(Spc.vgm, vgm("Gau","Exp"))
Spc.fit
SGau <- vgm(psill = 1727, model = "Gau", nugget = 443.97, range = 1.3998)
plot(Spc.vgm, model = SGau, col = "black", pch = 19, main = "SpC Variogram",
     ylab = "gamma", xlab = "distance (km)")

# Conductivity variogram
NO3dat <- dat[-which(is.na(dat$NO3.N.mgL)),]
NO3dat <- NO3dat[which(NO3dat$x<wwtp/1000),]
NO3.vgm <- variogram(NO3.N.mgL~1, loc = ~x+y, data = NO3dat, width = .2)
plot(NO3.vgm)

NH4.vgm <- variogram(NH4.N.mgL~1, loc = ~x+y, data = NO3dat, width = .2)
plot(NH4.vgm)

# cloud variogram
DO.vgm.cloud <- variogram(DO_mgL~1, loc= ~x+y, data=Spatialdat, cloud=T)
plot(DO.vgm.cloud, pch = 20, col = alpha("black", 0.2))
