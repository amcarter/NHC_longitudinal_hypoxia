###############################
# Process WQP data
# ACarter 4/27/2020

library(tidyverse)
library(streamMetabolizer)

####################################################################
# Load completed dataframe - sent it to Mike's computer to run
dat_piedmont <- readRDS("piedmontWQP_airP_done.rds")
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

plot(dat_piedmont$DOY, dat_piedmont$DO_persat*100, pch=21,
     col = "transparent",bg=alpha("black",.1),
     xlab = "Day of Year", ylab = "DO (% sat)", ylim = c(0,150), )
points(dat_piedmont$DOY[dat_piedmont$month=="October"], 100*dat_piedmont$DO_persat[dat_piedmont$month=="October"],
       col = alpha("brown3",.2), pch=20, bg="red")
boxplot(DO_persat*100 ~month, data=dat_piedmont, ylim = c(0,150), outline=FALSE,
        xlab = "Month", ylab = "DO (% sat)")
mtext("n=145175",side=1, line=-1.6, adj=.9)
mtext("1966-2020", side=1, line=-1.6, adj=.1)

w<- which(dat_piedmont$year%in%c(1966,1967,1968,2020))
tmp <- dat_piedmont[-w,]
boxplot(DO_persat*100 ~year, data=tmp, ylim = c(0,150), outline=FALSE,
        xlab = "Year", ylab = "DO (% sat)")
mtext("n=145175",side=1, line=-1.6, adj=.9)
mtext("1969-2019", side=1, line=-1.6, adj=.1)


abline(h=50, col="grey50")
w<- which(dat_piedmont$DO_persat<=.5)

length(w)/length(which(!is.na(dat_piedmont$DO_persat)))
seasonal <-dat_piedmont %>% group_by(month)%>%
  summarize(DO.persat.mean=mean(DO_persat, na.rm=T), 
            DO.mgl.mean = mean(DO_mgl, na.rm=T),
            n_hypox =length(which(DO_persat < .5))/length(month)*100)

par(mar = c(3,4,1,4))
boxplot(DO_persat*100 ~month, data=dat_piedmont, ylim = c(0,150), outline=FALSE,
        xlab = "", ylab = "")
mtext("DO (% sat)", side=2, line=-2, outer=T, adj=.55)
mtext("n=145175",side=3, line=-1.6, adj=.9)
mtext("1966-2020", side=3, line=-1.6, adj=.1)
par(new = T)
plot(seasonal$month, rep(-1000,12), ylim = c(0,25), axes=F, xlab="", ylab="")
polygon(c(seasonal$month, rev(seasonal$month)), c(seasonal$n_hypox, rep(0,12)), col = alpha("black", .2), border=F)
axis(4)
mtext("Pr(DO %sat < 50)", side=4, line=-2, outer=T, adj=.55)

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

