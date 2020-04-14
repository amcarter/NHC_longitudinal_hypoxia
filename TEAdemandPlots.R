################################
# Created by A Carter
# 10/31/2019

# Script for calculating metrics about DO time series 
# Plotting site characteristics

library(tidyverse)
library(lubridate)

#load element and flux data from sites
dat <- read.csv("data/FieldConcFluxes_NHCMC_20180702.csv", header=T, stringsAsFactors = F)
dat$datetime <- mdy_hm(dat$datetime)
glimpse(dat)

# remove the MC3_upper site and rename lower as MC3, rearrange based on site
dat <- dat[-8,]
dat$site[8] <- "MC3"
NHCsites2018 <- c('Mud','MC751','MC3','MC2','MC1','UNHC',
                  'NHC','NHC5','NHC4','NHC3','NHC2','NHC1')
dat$site <- factor(dat$site, levels = NHCsites2018)
dat <-dat[order(dat$site),]

# Remove high element concentrations below wwtp for correlations:
dat[dat$site=="NHC1",9:12]<- NA

#load DOmetric data:
DOmetrics <- read.csv("data/DOtimeseriesMetrics.csv", header=T, stringsAsFactors = F)

#Subset based on potential predictors:
DOmetPreds <- DOmetrics[,c(1,2,3,8,11,16:19)]
DOmetPreds <- cbind(DOmetPreds,dat[,c(2,5,6,7)])

plotvars <- function(DOmetPreds, variable){
 for(i in 2:ncol(DOmetPreds)){
         plot(DOmetPreds[,i], variable, ylab = F, xlab = colnames(DOmetPreds)[i])
         model <- lm(variable~DOmetPreds[,i])
         par(new=T)
         abline(model$coef[1], model$coef[2], lty=2)
         text(min(DOmetPreds[,i],na.rm=T),min(variable, na.rm=T), adj=c(0,0),
              labels = round(summary(model)$adj.r.squared,2))
         }       
}

par(mfrow = c(3,4), mar = c(4,2,0,0))
plotvars(DOmetPreds, dat$mean.DO_mgL.7.02.18)



# Build a data table with predictors and regression fits:

        
# Build a three panel figure with Storm O2 demand as x axis for each one
par(mfrow = c(3,1), mar = c(0,4,0,5), oma = c(4,0,1,0), xpd=T)

# Pull out the data for GHG fluxes
tmp <- dat[c(-8,-12,-13),]
tmp <- tmp[order(tmp$StormBOD.mgO2Ld),]

GHG <- rbind(tmp$flux.CO2.umolm2d,tmp$flux.CH4.umolm2d, tmp$flux.N2O.umolm2d)
rownames(GHG)<- c("CO2", "CH4","N2O")
colnames(GHG)<- tmp$StormBOD.mgO2Ld

TEAs <- rbind(tmp$DO.mg.L, tmp$NO3.N.mgL)
rownames(TEAs)<- c("O2", "NO3")
colnames(TEAs)<- tmp$StormBOD.mgO2Ld

EDs <- rbind(tmp$CO2.ugL, tmp$CH4.ugL, tmp$Fe.mgL)
 rownames(EDs)<- c("CO2", "CH4","Fe2+")
colnames(EDs)<- tmp$StormBOD.mgO2Ld

# Use this to convert to GHG equivalents
# GHG[2,]<-25*GHG[2,]
# GHG[3,]<- 298*GHG[3,]

barplot(GHG, border="white", space = 0.04, ylab = "flux (umol/m2/d)",
        col=c("grey20","grey50","red"), xaxt = "n")
legend("topright", inset=c(-.1,0), legend=c("CO2","CH4","N2O"), fill=c("grey20","grey50","red"), border="white", bty="n", cex=1.2)

plot(dat$StormBOD.mgO2Ld,dat$DO.mg.L, pch = 21, ylab = "Conc (mg/L)",  bg = "black", xaxt="n", cex=1.7)
par(new=T)
plot(dat$StormBOD.mgO2Ld[2:10], dat$NO3.N.mgL[2:10], pch = 21, bg="white", cex=1.7, xaxt = "n")

#points(dat$StormBOD.mgO2Ld[2:10],dat$SO4.mgl[2:10], pch = 19, col = "purple")
legend("topright", inset=c(-.12,0), legend = c("O2", "NO3"), pch =c(19,21), cex=1.2, bty = "n")

plot(dat$StormBOD.mgO2Ld,dat$CO2.ugL/1000, pch = 21, ylab = "Conc (mg/L)", bg="black", 
     xaxt="n", cex=1.7, xlab = "Ecosystem Oxygen Demand (mg/L/d)")
par(new=T)
plot(dat$StormBOD.mgO2Ld, dat$CH4.ugL, pch = 21,bg="white", cex=1.7)
par(new = T)
plot(dat$StormBOD.mgO2Ld, dat$Fe.mgL, pch = 21, col = "red", bg="red", cex=1.7)

legend("topright", inset=c(-.12,0), legend = c("CO2", "CH4", "Fe2+"), col = c("black","black", "red"), pch =c(19,21,19), cex=1.2, bty = "n")


# Remake plots with xaxis as BOM
# Build a three panel figure with Storm O2 demand as x axis for each one
par(mfrow = c(3,1), mar = c(0,4,0,5), oma = c(4,0,1,0), xpd=T)

# Pull out the data for GHG fluxes
tmp <- dat[c(-8,-11),]
tmp <- tmp[order(tmp$SedOrgM.per),]

GHG <- rbind(tmp$flux.CO2.umolm2d,tmp$flux.CH4.umolm2d, tmp$flux.N2O.umolm2d)
rownames(GHG)<- c("CO2", "CH4","N2O")
colnames(GHG)<- tmp$SedOrgM.per


barplot(GHG, border="white", space = 0.04, ylab = "flux (umol/m2/d)",
        col=c("grey20","grey50","red"), xaxt = "n")
legend("topright", inset=c(-.1,0), legend=c("CO2","CH4","N2O"), fill=c("grey20","grey50","red"), border="white", bty="n", cex=1.2)

plot(dat$SedOrgM.per,dat$DO.mg.L, pch = 21, ylab = "Conc (mg/L)",  bg = "black", xaxt="n", cex=1.7)
par(new=T)
plot(dat$SedOrgM.per[2:10], dat$NO3.N.mgL[2:10], pch = 21, bg="white", cex=1.7, xaxt = "n")

#points(dat$StormBOD.mgO2Ld[2:10],dat$SO4.mgl[2:10], pch = 19, col = "purple")
legend("topright", inset=c(-.12,0), legend = c("O2", "NO3"), pch =c(19,21), cex=1.2, bty = "n")

plot(dat$SedOrgM.per,dat$CO2.ugL/1000, pch = 21, ylab = "Conc (mg/L)", bg="black", 
     xaxt="n", cex=1.7, xlab = "Ecosystem Oxygen Demand (mg/L/d)")
par(new=T)
plot(dat$SedOrgM.per, dat$CH4.ugL, pch = 21,bg="white", cex=1.7)
par(new = T)
plot(dat$SedOrgM.per, dat$Fe.mgL, pch = 21, col = "red", bg="red", cex=1.7)

legend("topright", inset=c(-.12,0), legend = c("CO2", "CH4", "Fe2+"), col = c("black","black", "red"), pch =c(19,21,19), cex=1.2, bty = "n")


