#####################
# Conceptual figure for NHC hypoxia paper
# Across site variation:
# Panel 1: single pre-storm day of data from all 12 sites
#       2: Storm recovery across sites
#       3: Seasonal variation across sites
# Within site variation
# Panel 4: all days of data for NHC 2018
#       5: storm recovery trajectories for NHC
#       6: Interannual variability at one site
#   
library(lubridate)
library(zoo)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

################################################
# Panel 1: single pre-storm day of data from all 12 sites

# Read in DO data and daily site summarys:
# DOSiteSummary <- read.csv("data/DailyDOsummary.csv", header=T)
DOdat <- read.csv("data/raw/allSites_20180706.csv", header = T, stringsAsFactors = F)
NHCsites2018 <- c('Mud','MC751','MC3','MC2','MC1','UNHC',
                  'NHC','NHC5','NHC4','NHC3','NHC2','NHC1')
DOdat$site <- factor(DOdat$site, levels = NHCsites2018)
DOdat$DateTime_UTC <- ymd_hms(DOdat$DateTime_UTC)
DOdat$DateTime <- with_tz(DOdat$DateTime_UTC, tz="EST")
# The targeted storm occured on Jun 26:
DOdat <- select(DOdat, persatDO, DateTime, site)
DOdat$persatDO <- na.approx(DOdat$persatDO)
DOdat$persatDO.smooth<- rollmean(DOdat$persatDO, 24,na.pad=T)



startdate <- ymd_hms("2018-06-25 00:00:00", tz="EST")
enddate <- ymd_hms("2018-06-26 00:00:00", tz="EST")

pal <- colorRampPalette(c("brown3",
                          "darkblue"))
cols <- pal(12)

site.deets <- read_csv("2018 NHC MC Longitudinal/Longitudinal Survey May 2018/Longitudinal Data/NHC_sensorsite_data.csv")

png("figures/Figure5_six_panel_data_fig_supplement.png", width = 8, height = 5,
    units = "in", res = 300)
#################################################
# Parameters for multipanel graph

v<-c(rep(c(rep(1,15),2,rep(3,22),4,rep(5,15)),15),
     rep(6,54),
     rep(c(rep(7,15),8,rep(9,22),10,rep(11,15)),15))
m <- matrix(v, ncol=54, byrow=T)

layout(m)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 2.2, 2))

# text and line colors and thicknesses
bor.col <- "grey60"
text.col <- "grey20"
lab.col <- "black"

plot(1,xlim = c(startdate,enddate), ylim = c(0, 120), axes=F,
     type = "n")
mtext(paste(letters[1], "   12 sites", sep=" "), side = 3, 
      line = -1.5, #how far from the side should it be, 
      adj = 0.02, # percent of way along the side to print
      cex = 1, col = text.col)
axis(2, col = bor.col, col.axis = text.col, at = seq(0, 100, 25))
box(col = bor.col)
# plot all
for(i in 1:length(NHCsites2018)){
  lines(DOdat[DOdat$site==NHCsites2018[i],]$DateTime,
        100*DOdat[DOdat$site==NHCsites2018[i],]$persatDO.smooth,
        col = cols[i], lwd = 2)
}

plot(1,axes=F, type="n")
####################################################################
# Panel 2: Storm recovery across sites

startdate <- ymd_hms("2018-06-25 0:00:00", tz="EST")
enddate <- ymd_hms("2018-07-02 24:00:00", tz="EST")
DOdat$persatDO.smooth<- rollmean(DOdat$persatDO, 36, na.pad=T)
# Plot up:
plot(1,xlim = c(startdate,enddate), ylim = c(0, 120), axes=F, type="n")
mtext(paste(letters[2], "   12 sites", sep=" "), side = 3, line = -1.5, #how far from the side should it be, 
      adj = 0.02, # percent of way along the side to print
      cex = 1, col = text.col)
box(col = bor.col)
#plot high
for(i in 1:length(NHCsites2018)){
  lines(DOdat[DOdat$site==NHCsites2018[i],]$DateTime, 
        100*DOdat[DOdat$site==NHCsites2018[i],]$persatDO.smooth,
        col = cols[i], lwd = 2)
}
plot(1,axes=F, type="n")

# plot(DOdat[DOdat$site=="MC3",]$DateTime, 100*DOdat[DOdat$site=="MC3",]$persatDO.smooth,
#   xlim = c(startdate,enddate), ylim = c(0, 120), type = "l")
# for(i in 1:length(NHCsites2018)){
#   lines(DOdat[DOdat$site==NHCsites2018[i],]$DateTime, 
#         100*DOdat[DOdat$site==NHCsites2018[i],]$persatDO.smooth,
#         col = cols[i], lwd = 2)
# }

#############################################################
# Panel 3: Seasonal variation across sites

dat <- read.csv("data/raw/2019SPsites.csv", header = T, stringsAsFactors = F)
dat$DateTime_UTC <- ymd_hms(dat$DateTime_UTC)
dat$DateTime <- with_tz(dat$DateTime_UTC, tz="EST")
dat$Date <- as.Date(dat$DateTime, tz = "EST")


# Summarize by day
datdaily <- dat %>% dplyr::group_by(site, Date) %>% 
  dplyr::summarise(mean.persatDO = mean(persatDO, na.rm=T),
                   sd.persatDO = sd(persatDO, na.rm=T)) 
datdaily$mean.persatDO <- na.approx(datdaily$mean.persatDO, na.rm=F)

# Calculate exceedence curve for a site

plot_exceedence <- function(tmp, col="black", lty=1){
  tmp <- tmp[which(!is.na(tmp$mean.persatDO)),]
  tmp$mean.persatDO <- sort(tmp$mean.persatDO, decreasing=T)
  tmp$index <- seq(1:nrow(tmp))
  tmp$freq <- 100*(tmp$index/(1+nrow(tmp)))
  lines(tmp$freq,tmp$mean.persatDO*100, lwd = 2, col = col)
}

NHCsites <- c(1,2,6,7)

plot(1,xlim = c(0,100), ylim = c(0,120), type = "n",axes=F)
mtext(paste(letters[3], "   4 sites", sep=" "), side = 3, line = -1.5, #how far from the side should it be, 
      adj = 0.02, # percent of way along the side to print
      cex = 1, col = text.col)
box(col = bor.col)


for(i in NHCsites){
  plot_exceedence(datdaily[datdaily$site==NHCsites2018[i], "mean.persatDO"], 
                  cols[i])
}

plot(1,axes=F, type="n")

#########################################################
# Within site variation
# Panel 4: all days of data for NHC 2018

dat <- read_csv(file = "data/raw/2019SPsites.csv")
dat$DateTime <- with_tz(dat$DateTime_UTC, tz = "EST")
dat$Date <- as.Date(dat$DateTime, tz = "EST")

NHC <- dat[dat$site=="NHC",]%>% select(-site, -DateTime_UTC)
NHC$persatDO.smooth <- rollmean(NHC$persatDO, 32, na.pad=T)
NHC <- NHC %>% 
  mutate(year = year(Date)) %>%
  filter(year == 2019)
# NHC$col[NHC$month==3 |NHC$month==2]<- "darkblue"

plot(NHC[NHC$Date==as.Date("2019-02-01"),]$DateTime,
     100*NHC[NHC$Date==as.Date("2019-02-01"),]$persatDO,
     ylim = c(0,120), type = "n", col.axis = text.col,
     yaxt="n", xlab = "Time")
mtext(paste(letters[4], "   365 days", sep=" "), side = 3, line = -1.5, #how far from the side should it be, 
      adj = 0.02, # percent of way along the side to print
      cex = 1, col = text.col)
axis(2, col =bor.col, col.axis = text.col, at = seq(0, 100, 25))
mtext("time", side = 1, outer = TRUE, cex = 1, line = 2.2, adj = .15, col = lab.col)


for(i in unique(NHC$Date)){
  par(new=T)
  plot(NHC[NHC$Date==i,]$DateTime, 100*NHC[NHC$Date==i,]$persatDO,
       type="l",col = alpha("black", 0.2), lwd = 2, 
       ylab="",xlab="",xaxt="n", yaxt="n",ylim = c(0,120))
}
box(col = bor.col)

plot(1,axes=F, type="n")

# Panel 5: storm recovery trajectories for NHC
# library(dygraphs)
# library(xts)
# NHC.dy <- xts(x=select(NHC, persatDO, Level_m), order.by=NHC$DateTime)
#   
# dygraph(NHC.dy)%>%dyRangeSelector()


stormdates <- as.Date(c("2019-01-13",#"2019-01-20","2019-01-24","2019-02-13",
                        "2019-02-16",#"2019-02-18","2019-02-21","2019-02-24",
                        "2019-03-21","2019-04-06",#"2019-04-09","2019-04-13",
                        "2019-04-20","2019-06-08",#"2019-06-19",
                        "2019-07-23",
                        "2019-08-05","2019-08-14","2019-10-20","2019-11-24",
                        "2019-12-01"))#,"2019-12-14", "2019-12-24"))

# all 
# stormdates <- as.Date(c("2019-01-13","2019-01-20","2019-01-24","2019-02-13",
#                         "2019-02-16","2019-02-18","2019-02-21","2019-02-24",
#                         "2019-03-21","2019-04-06","2019-04-09","2019-04-13",
#                         "2019-04-20","2019-06-08","2019-06-19","2019-07-23",
#                         "2019-08-05","2019-08-14","2019-10-20","2019-11-24",
#                         "2019-12-01","2019-12-14", "2019-12-24"))

stormdates <- stormdates - 2
endinterval <- stormdates+9

plot(1, type = "n", ylim = c(0,120), xlim = c(-2.5,6.5), axes=F)
mtext(paste(letters[5], "    12 events", sep=" "), side = 3, line = -1.5, #how far from the side should it be, 
      adj = 0.02, # percent of way along the side to print
      cex = 1, col = text.col)
axis(1, col = bor.col, col.axis =text.col, at = seq(-2, 6, 1) )
box(col = bor.col)
mtext("days post storm", side = 1, outer = TRUE, cex = 1, line = 2.2, col = lab.col)



for(i in 1:length(stormdates)){
  par(new = T)
  plot(NHC[NHC$Date >= stormdates[i] & NHC$Date<endinterval[i], ]$DateTime,
       100*NHC[NHC$Date >=stormdates[i]&NHC$Date<endinterval[i],]$persatDO.smooth,
       ylim = c(0,120), type = "l",
       col = alpha(1, 0.5),lwd = 2,
       xaxt="n",yaxt="n", xlab = '', ylab = '')
}
box(col = bor.col)

plot(1,axes=F, type="n")

##########################################################
# Panel 6: Interannual variability at one site

dat <- read_csv("data/raw/NHCdat.csv")
dat$DateTime <- with_tz(dat$DateTime_UTC, tz="EST")
dat$Date <- as.Date(dat$DateTime, tz = "EST")

datdaily <- dat %>% dplyr::group_by(Date) %>% 
  dplyr::summarise(mean.persatDO = mean(persatDO, na.rm=T),
                   sd.persatDO = sd(persatDO, na.rm=T)) 
datdaily$mean.persatDO <- na.approx(datdaily$mean.persatDO)
datdaily$year <- strftime(datdaily$Date, format = "%Y")


plot(1,xlim = c(0,100), ylim = c(0,120), type = "n",axes=F)
mtext(paste(letters[6], "   3 years", sep=" "), side = 3, line = -1.5, #how far from the side should it be, 
      adj = 0.02, # percent of way along the side to print
      cex = 1, col = text.col)
axis(1, col = bor.col, col.axis = text.col, at = seq(0,100, 20) )
box(col = bor.col)
#legend("bottomleft", legend = c("2017","2018","2019"),cex=.7,bty="n", lty = c(1,2,3),col=alpha("black", .7), lwd=2)

for(i in c(2017,2018,2019)){
  plot_exceedence(datdaily[datdaily$year==i, "mean.persatDO"], col=alpha("black", .7),lty=i)
}

mtext("% exceedence", side = 1, outer = TRUE, cex = 1, line = 2.2, adj = .9, col = lab.col)
mtext("DO (% saturation)", side = 2, outer = TRUE, cex = 1, line = 2.2, col = lab.col)
mtext("Within Site",side = 4, outer=TRUE, cex=1, line = .3,adj = .2, col = lab.col )
mtext("Between Sites",side = 4, outer=TRUE, cex=1, line = .3,adj=.8, col = lab.col )
par(mfrow = c(1,1), oma = c(0,0,0,0), new = T)
plot(1,axes = F, type = "n")
legend("top", legend = c("Upstream Site", "Downstream Site"),
       col = c(cols[1],cols[12]),
       lty = c(1,1), ncol = 2, bty = "n", lwd = 2,
       cex = .8)
dev.off()
