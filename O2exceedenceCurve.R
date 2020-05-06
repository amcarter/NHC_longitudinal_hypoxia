#########################################
# Dissolved Oxygen exceedence curves - the sagginess index

# A Carter 3/10/2020

# Calculate and plot flow exceedence curves for annual DO data

library(lubridate)

dat <- read.csv("data/raw/2019SPsites.csv", header = T, stringsAsFactors = F)
dat$DateTime_UTC <- ymd_hms(dat$DateTime_UTC)
dat$DateTime <- with_tz(dat$DateTime_UTC, tz="EST")
dat$Date <- as.Date(dat$DateTime, tz = "EST")


# Summarize by day
datdaily <- dat %>% dplyr::group_by(site, Date) %>% 
  dplyr::summarise(mean.persatDO = mean(persatDO, na.rm=T),
                   sd.persatDO = sd(persatDO, na.rm=T)) 
datdaily$mean.persatDO <- na.approx(datdaily$mean.persatDO)

# Calculate exceedence curve for a site

plot_exceedence <- function(tmp, col="black", lty=1){
  tmp <- tmp[which(!is.na(tmp))]
  tmp<- sort(tmp, decreasing=T)
  index <- seq(1:length(tmp))
  freq <- 100*(index/(1+length(tmp)))
  lines(freq,tmp, lwd = 2, col = col, lty=lty)
}

NHCsites <- c("MC751", "Mud", "UNHC", "NHC")
cols <- c("brown3", colorRampPalette(c("grey20", "grey80"))(3))

par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(datdaily$mean.persatDO, datdaily$mean.persatDO, 
     xlim = c(0,100), ylim = c(0,120), type = "n",
     ylab = "Daily mean DO (% sat)", xlab = "Exceedence Frequency", 
     main = "2019")

for(i in 1:4){
  plot_exceedence(100*datdaily[datdaily$site==NHCsites[i],]$mean.persatDO, cols[i])
}

###################################
# Plot curves for years in NHC

dat <- read.csv("data/raw/NHCdat.csv", header = T, stringsAsFactors = F)
dat$DateTime_UTC <- ymd_hms(dat$DateTime_UTC)
dat$DateTime <- with_tz(dat$DateTime_UTC, tz="EST")
dat$Date <- as.Date(dat$DateTime, tz = "EST")

datdaily <- dat %>% dplyr::group_by(Date) %>% 
  dplyr::summarise(mean.persatDO = mean(persatDO, na.rm=T),
                   sd.persatDO = sd(persatDO, na.rm=T),
                   Q=mean(discharge_cms, na.rm=T)) 
datdaily$mean.persatDO <- na.approx(datdaily$mean.persatDO)
plot(datdaily$Date, datdaily$Q, log="y", type="l")

datdaily$year <- strftime(datdaily$Date, format = "%Y")


par(mfrow = c(1,2), mar = c(4,4,1,1))
plot(datdaily$mean.persatDO, datdaily$mean.persatDO, 
     xlim = c(0,100), ylim = c(0,120), type = "n",
     ylab = "Daily mean DO (% sat)", xlab = "Exceedence Frequency", 
     main = "")
legend("bottomleft", legend = c("2017","2018","2019"),cex=.7,bty="n", lty = c(1,2,3), lwd=2)
for(i in c(2017,2018,2019)){
  plot_exceedence(100*datdaily[datdaily$year==i,]$mean.persatDO, lty=i)
}



###################################
# DO rating curve

plot(datdaily$mean.persatDO, datdaily$Q, log="y")
datdaily$Q <- na.approx(datdaily$Q)
par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(datdaily$Q, datdaily$mean.persatDO, 
     xlim = c(0,100), ylim = c(-10,10), type = "n",
     ylab = "Q", xlab = "Exceedence Frequency", 
     main = "")

for(i in c(2017,2018,2019)){
  plot_exceedence(log(datdaily[datdaily$year==i,]$Q), lty=i)
}
for(i in c(2017,2018,2019)){
  par(new=T)
  plot(datdaily[datdaily$year==i,]$Date, log(datdaily[datdaily$year==i,]$Q), ylim = c(-10,10), axes=FALSE,
       xlab="", ylab = "", type="l", col=i)
}

plot_exceedence(log(datdaily$Q))
