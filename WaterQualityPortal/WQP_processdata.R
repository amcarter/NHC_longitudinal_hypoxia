###############################
# Process WQP data
# ACarter 4/27/2020

library(tidyverse)
library(streamMetabolizer)
setwd(hypox_projdir)
####################################################################
# Load completed dataframe - sent it to Mike's computer to run
dat_piedmont <- read_csv("WaterQualityPortal/WQP_NHDlinked.csv")
length(which(is.na(dat_piedmont$airPres_kPa)))/nrow(dat_piedmont)
plot(dat_piedmont$temp_C, dat_piedmont$air_temp_C)


###########################################################
# Convert DO data to percent saturation using StreamMetabolizer package
# Pressure must be converted to millibars first: 1 kPa = 10 mbar
dat_piedmont$DO_sat<- calc_DO_sat(dat_piedmont$temp_C, 10*dat_piedmont$airPres_kPa)
dat_piedmont$DO_persat <- dat_piedmont$DO_mgl/dat_piedmont$DO_sat 
dat_piedmont$month <- factor(strftime(dat_piedmont$DateTime_UTC, "%b"), 
                             levels = format(ISOdatetime(2000,1:12,1,0,0,0),"%b"))
dat_piedmont$DOY <- strftime(dat_piedmont$DateTime_UTC, format = "%j")

summary(dat_piedmont)

plot(dat_piedmont$DOY, dat_piedmont$DO_persat*100, pch=21,
     col = "transparent",bg=alpha("black",.1),
     xlab = "Day of Year", ylab = "DO (% sat)", ylim = c(0,150), )
points(dat_piedmont$DOY[dat_piedmont$month=="October"], 100*dat_piedmont$DO_persat[dat_piedmont$month=="October"],
       col = alpha("brown3",.2), pch=20, bg="red")
boxplot(DO_persat*100 ~month, data=dat_piedmont, ylim = c(0,150), outline=FALSE,
        xlab = "Month", ylab = "DO (% sat)")
mtext("n=145175",side=1, line=-1.6, adj=.9)
mtext("1966-2020", side=1, line=-1.6, adj=.1)

################################################
# final plot
w<- which(dat_piedmont$DO_persat<=.5)
which(is.na(dat_piedmont$month))

length(w)/length(which(!is.na(dat_piedmont$DO_persat)))
seasonal <-dat_piedmont %>% group_by(month)%>%
  summarize(DO.persat.mean=mean(DO_persat, na.rm=T), 
            DO.mgl.mean = mean(DO_mgl, na.rm=T),
            n_hypox =length(which(DO_persat < .5))/length(month)*100)

slope.75<-quantile(dat_piedmont$slope, .75, na.rm=T)
seasonal.highslope <-dat_piedmont[dat_piedmont$slope>slope.75,] %>% group_by(month)%>%
  summarize(DO.persat.mean=mean(DO_persat, na.rm=T), 
            DO.mgl.mean = mean(DO_mgl, na.rm=T),
            n_hypox =length(which(DO_persat < .5))/length(month)*100)

slope.25 <- quantile(dat_piedmont$slope, .25, na.rm=T)
seasonal.lowslope <-dat_piedmont[dat_piedmont$slope<slope.25,] %>% group_by(month)%>%
  summarize(DO.persat.mean=mean(DO_persat, na.rm=T), 
            DO.mgl.mean = mean(DO_mgl, na.rm=T),
            n_hypox =length(which(DO_persat < .5))/length(month)*100)


png("figures/WQP_DO_seasonal.png", width=5, height=4, units="in", res=300)

par(mar = c(3,4,1,4), oma = c(0,0,3,0))

boxplot(DO_persat*100 ~month, data=dat_piedmont, ylim = c(0,150), outline=FALSE,
        xlab = "", ylab = "", cex.axis=.8)
mtext("DO (% sat)", side=2, line=-2, outer=T, adj=.55)
mtext("n=145175",side=3, line=-1.6, adj=.9)
mtext("1966-2020", side=3, line=-1.6, adj=.1)
par(new = T)


plot(seasonal$month, rep(-1000,12), ylim = c(0,25), axes=F, xlab="", ylab="")
polygon(c(seq(1,12),rev(seq(1,12))), c(seasonal.lowslope$n_hypox[1:12], rep(0,12)), col = alpha("black", .2), border=F)
polygon(c(seq(1,12),rev(seq(1,12))), c(seasonal.highslope$n_hypox[1:12], rep(0,12)), col = alpha("black", .4), border=F)
axis(4, cex.axis=.8)
mtext("Pr(DO % sat < 50)", side=4, line=-2, outer=T, adj=.55)
par(new=T, oma=c(0,0,0,0))

legend("top", legend=c(paste("Lower 25% of slope\n     (<",round(slope.25*1000, 1),"m/km)"),
                       paste("Upper 25% of slope\n     (>",round(slope.75*1000, 1),"m/km)")), 
       cex=.8, bty="n",ncol=2, xpd=NA,
       fill=c(alpha("black",.2),alpha("black",.6)),
       border=c(NA,NA))

dev.off()
#############################################################
#quantile regression on lower 95% of slope data
maxslope <- quantile(dat_piedmont$slope, .95, na.rm=T)

dat_piedmont$hypox <- 0
dat_piedmont$hypox[dat_piedmont$DO_persat<.5]<-1
dat_piedmont$hypox[is.na(dat_piedmont$DO_persat)]<-NA

plot(dat_piedmont$slope, dat_piedmont$DO_persat, xlim = c(0,maxslope))
tmp <- dat_piedmont[dat_piedmont$slope<maxslope,]
tmp$slope <- tmp$slope*1000
tmp$DO_persat <- tmp$DO_persat*100
mm<- glm(DO_persat~slope,data=tmp )

confint(mm)
summary(mm)
abline(.7973, 12.66, col=2)

dat_piedmont_ordered <- dat_piedmont[order(dat_piedmont$slope),]
dat_piedmont_ordered <- dat_piedmont_ordered[which(!is.na(dat_piedmont_ordered$slope)),]
dat_piedmont_ordered <- dat_piedmont_ordered[which(!is.na(dat_piedmont_ordered$hypox)),]
dat_piedmont_ordered <- dat_piedmont_ordered[dat_piedmont_ordered$month %in% c("Jul","Aug","Sep"),]

l <- nrow(dat_piedmont_ordered)
nbins <- 10
bins<-rep(1:nbins, each=l/nbins)
bins <- c(bins,rep(10,8))
avg_hpx <- data.frame(slope=tapply(dat_piedmont_ordered$slope,bins, mean),
                      hpx=tapply(dat_piedmont_ordered$hypox,bins, mean),
                      hpx.sd=tapply(dat_piedmont_ordered$hypox,bins, sd),
                      hpx.975=tapply(dat_piedmont_ordered$hypox,bins, quantile, .975),
                      DO=tapply(dat_piedmont_ordered$DO_persat,bins, mean))
plot(avg_hpx$slope, avg_hpx$DO, ylim = c(0,1.2), type="b", lwd=2, xlim = c(0,0.01))                      
polygon(c(avg_hpx$slope, rev(avg_hpx$slope)), c(avg_hpx$hpx+avg_hpx$hpx.sd, rev(avg_hpx$hpx-avg_hpx$hpx.sd)), 
        col=alpha("black", .3), border=NA)
                      

xrng = c(0,maxslope)
cutseq = seq(xrng[1], xrng[2], length.out=11)
bins = cut(tmp$slope, cutseq)
out<- data.frame(slope=cutseq[1:10],
                 hpx=tapply(tmp$hypox, bins, mean, na.rm=T),
                 hpx.sd=tapply(tmp$hypox, bins, sd, na.rm=T))

plot(out$slope, out$hpx, ylim = c(0,.5),type="l", lwd=2)
polygon(c(out$slope, rev(out$slope)), c(rep(0,10),rev(out$hpx+out$hpx.sd)), col=alpha("black", .3), border=NA)

#################################################################

w<- which(dat_piedmont$year%in%c(1966,1967,1968,2020))
dat_piedmont <- dat_piedmont[-w,]
boxplot(DO_persat*100 ~year, data=tmp, ylim = c(0,150), outline=FALSE,
        xlab = "Year", ylab = "DO (% sat)")
mtext("n=145175",side=1, line=-1.6, adj=.9)
mtext("1969-2019", side=1, line=-1.6, adj=.1)



abline(h=50, col="grey50")


###################################################################################################
# Change in DO distribution over time
tmp <- dat_piedmont%>% group_by(year, month)%>%
  summarise(DO_persat_mean =mean(DO_persat, na.rm=T))
plot(tmp$year[tmp$month=="Jul"], tmp$DO_persat_mean[tmp$month=="Jul"], type="l",lty=1)
lines(tmp$year[tmp$month=="Aug"], tmp$DO_persat_mean[tmp$month=="Aug"], lty=2)
lines(tmp$year[tmp$month=="Oct"], tmp$DO_persat_mean[tmp$month=="Oct"], lty=3)
plot(tmp$year[tmp$month=="Jul"], tmp$DO_persat_mean[tmp$month=="Jul"]/tmp$DO_persat_mean[tmp$month=="Oct"], type="l")


# Violin Plots
library(vioplot)
x1 <- dat_piedmont$DO_persat[dat_piedmont$month==levels(dat_piedmont$month)[1]]
x2 <- dat_piedmont$DO_persat[dat_piedmont$month==levels(dat_piedmont$month)[2]]
x3 <- dat_piedmont$DO_persat[dat_piedmont$month==levels(dat_piedmont$month)[3]]
x4 <- dat_piedmont$DO_persat[dat_piedmont$month==levels(dat_piedmont$month)[4]]
x5 <- dat_piedmont$DO_persat[dat_piedmont$month==levels(dat_piedmont$month)[5]]
x6 <- dat_piedmont$DO_persat[dat_piedmont$month==levels(dat_piedmont$month)[6]]
x7 <- dat_piedmont$DO_persat[dat_piedmont$month==levels(dat_piedmont$month)[7]]
x8 <- dat_piedmont$DO_persat[dat_piedmont$month==levels(dat_piedmont$month)[8]]
x9 <- dat_piedmont$DO_persat[dat_piedmont$month==levels(dat_piedmont$month)[9]]
x10 <- dat_piedmont$DO_persat[dat_piedmont$month==levels(dat_piedmont$month)[10]]
x11 <- dat_piedmont$DO_persat[dat_piedmont$month==levels(dat_piedmont$month)[11]]
x12 <- dat_piedmont$DO_persat[dat_piedmont$month==levels(dat_piedmont$month)[12]]

vioplot(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12,
        names=levels(dat_piedmont$month),
        col="grey60")
title("Violin Plots of Miles Per Gallon")


plot(dat_piedmont$DO_mgl, dat_piedmont$DO_persat)

annual <- dat_piedmont %>%group_by(year)%>%
  summarize(DO_persat.mean=mean(DO_persat, na.rm=T),
            DO_persat.sd=sd(DO_persat, na.rm=T))
plot(annual$year, annual$DO_persat.mean, ylim = c(0,1.2), type="l")
polygon(x=na.approx(c(annual$year, rev(annual$year))), 
        y=na.approx(c(annual$DO_persat.mean+annual$DO_persat.sd, 
                      rev(annual$DO_persat.mean-annual$DO_persat.sd))),
        col=alpha("black", .2), bty="n")
abline(h=1)

hist(dat_piedmont$DO_persat, xlim = c(0,1.5), ylim = c(0,36000))
par(new=T)
hist(dat_piedmont$DO_persat[dat_piedmont$month=="October"], xlim = c(0,1.5), 
     ylim = c(0,36000), col = "brown3", axes=F)

plot(dat_piedmont$DateTime_UTC, dat_piedmont$DO_persat, col = alpha("black",0.2), pch=19, ylim = c(0,1.5))
points(dat_piedmont$DateTime_UTC[dat_piedmont$month=="October"],
       dat_piedmont$DO_persat[dat_piedmont$month=="October"],pch=19, col=alpha("brown3",.6) )
plot(annual$year, annual$DO_persat.mean, type="l", ylim = c(0,1.5))



ncSummary <- DOdat %>%
  group_by(MonitoringLocationIdentifier) %>%
  summarise(count=n(),
            max = max(value, na.rm = TRUE)) %>%
  filter(count > 10) %>%
  arrange(-count) %>%
  left_join(siteInfo, by = "MonitoringLocationIdentifier")


col_types <- c("darkblue","dodgerblue","green4","gold1","orange","brown","red")
leg_vals <- unique(as.numeric(quantile(ncSummary$max,
                                       probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,.99,1), na.rm=TRUE)))
pal = colorBin(col_types, ncSummary$max, bins = leg_vals)
rad <-3*seq(1,4,length.out = 16)
ncSummary$sizes <- rad[as.numeric(cut(ncSummary$count, breaks=16))]

leaflet(data=ncSummary) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(~dec_lon_va,~dec_lat_va,
                   fillColor = ~pal(max),
                   radius = ~sizes,
                   fillOpacity = 0.8, opacity = 0.8,stroke=FALSE,
                   popup=~station_nm) %>%
  addLegend(position = 'bottomleft',
            pal=pal,
            values=~max,
            opacity = 0.8,
            labFormat = labelFormat(digits = 1),
            title = 'Max Value')

