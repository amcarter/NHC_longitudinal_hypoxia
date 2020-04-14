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

plot_exceedence <- function(tmp, col){
  tmp <- tmp[which(!is.na(tmp$mean.persatDO)),]
  tmp$mean.persatDO <- sort(tmp$mean.persatDO, decreasing=T)
  tmp$index <- seq(1:nrow(tmp))
  tmp$freq <- 100*(tmp$index/(1+nrow(tmp)))
  lines(tmp$freq,tmp$mean.persatDO*100, lwd = 2, col = col)
}

NHCsites <- c("MC751", "Mud", "UNHC", "NHC")
cols <- c("brown3", colorRampPalette(c("grey20", "grey80"))(3))

par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(datdaily$mean.persatDO, datdaily$mean.persatDO, 
     xlim = c(0,100), ylim = c(0,120), type = "n",
     ylab = "Daily mean DO (% sat)", xlab = "Exceedence Frequency", 
     main = "2019")

for(i in 1:4){
  plot_exceedence(datdaily[datdaily$site==NHCsites[i], "mean.persatDO"], cols[i])
}

###################################
# Plot curves for years in NHC

dat <- read.csv("data/raw/NHCdat.csv", header = T, stringsAsFactors = F)
dat$DateTime_UTC <- ymd_hms(dat$DateTime_UTC)
dat$DateTime <- with_tz(dat$DateTime_UTC, tz="EST")
dat$Date <- as.Date(dat$DateTime, tz = "EST")

datdaily <- dat %>% dplyr::group_by(Date) %>% 
  dplyr::summarise(mean.persatDO = mean(persatDO, na.rm=T),
                   sd.persatDO = sd(persatDO, na.rm=T)) 
datdaily$mean.persatDO <- na.approx(datdaily$mean.persatDO)
datdaily$year <- strftime(datdaily$Date, format = "%Y")


par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(datdaily$mean.persatDO, datdaily$mean.persatDO, 
     xlim = c(0,100), ylim = c(0,120), type = "n",
     ylab = "Daily mean DO (% sat)", xlab = "Exceedence Frequency", 
     main = "NHC")

for(i in c(2017,2018,2019)){
  plot_exceedence(datdaily[datdaily$year==i, "mean.persatDO"], alpha("black", 0.5))
}
