### Plot average annual hydrograph for USGS Blands Station 
# 2020-02-19
### Plot average annual hydrograph for SP NHC station 
# 2020-10-26

# Updates:
#       2020-03-27 Added section to calculate NHC baseflow and compare to wwtp discharge
#       2020-04-06 Updated with baseflow separation package RHydro
library(zoo)
library(ggplot2)
#install.packages("RHydro", repos="http://R-Forge.R-project.org")
library(RHydro)
library(StreamPULSE)
setwd(hypox_projdir)
dat<- read.csv("data/raw/NHC_Blands_annualdata.csv", header=T)
dat <- read_csv("data/NHC_Q_and_temp.csv")
dat$date <- as.Date(dat$DateTime_UTC, tz = "EST")
dat <- dat %>% group_by(date) %>%
        summarise(discharge_m3s = mean(discharge_m3s, na.rm = T),
                  temp_C = mean(water_temp_C, na.rm = T)) 
dat$doy <- format(dat$date, "%j")
dat$year <- year(dat$date)
nhc2018 <- dat %>%
        filter(year == 2018)

nhc_sum <- dat %>%
        group_by(doy) %>%
        summarise(q.25 = quantile(discharge_m3s, 0.025, na.rm = T),
                  q.75 = quantile(discharge_m3s, 0.975, na.rm = T),
                  t.25 = quantile(temp_C, 0.025, na.rm = T),
                  t.75 = quantile(temp_C, 0.975, na.rm = T))

nhc_day <- nhc2018 %>% filter(date == as.Date("2018-05-25"))
nhc_sens <- nhc2018 %>% 
        filter(date >= as.Date("2018-06-12"), date <= as.Date("2018-07-06"))
summary(dat)
png("figures/NHC_average_Q_and_t.png", width = 7, height = 5, 
    res = 300, units = "in")
        par(mfrow = c(2,1), oma = c(0,0,2,0), mar = c(0,4,1,2))                  
        plot(nhc_sum$doy, nhc_sum$q.75, log = 'y', type = "n",
             ylab = "discharge (m3/s)", ylim = c(.0005,200), 
             xaxt = "n", xlab = "") 
        polygon(x = na.approx(c(nhc_sum$doy, rev(nhc_sum$doy))),
                y = na.approx(c(nhc_sum$q.25, rev(nhc_sum$q.75))),
                col = "grey80", border = NA)
        lines(nhc2018$doy, nhc2018$discharge_m3s, lwd = 2)
        lines(nhc_sens$doy, nhc_sens$discharge_m3s, lwd = 2, col = "brown3")
        
        points(nhc_day$doy, nhc_day$discharge_m3s, pch = 19, 
               col = "brown3", cex = 1.3)
        
        par(mar = c(4,4,0,2))
        plot(nhc_sum$doy, nhc_sum$t.75, type = "n",
             ylab = "water temp (C)", ylim = c(0,28), xlab = "day of year") 
        polygon(x = na.approx(c(nhc_sum$doy, rev(nhc_sum$doy))),
                y = na.approx(c(nhc_sum$t.25, rev(nhc_sum$t.75))),
                col = "grey80", border = NA)
        lines(nhc2018$doy, nhc2018$temp_C, lwd = 2)
        lines(nhc_sens$doy, nhc_sens$temp_C, lwd = 2, col = "brown3")
        points(nhc_day$doy, nhc_day$temp_C, pch = 19, col = "brown3", cex = 1.3)
        
        
        par(new = T, oma = c(0,0,0,0), mfrow = c(1,1))
        plot(1,1, type = "n", xaxt = "n", yaxt = "n",
             xlab = "", ylab = "", bty = "n" )
        legend("top", ncol = 2, bty = "n", cex = .9,
               legend = c("Longitudinal sampling","sensor deployment", 
                          "   2018 data", "2016-2019 IQ range"),
               lty = c(NA, 1, 1, NA),
               pch = c(19, NA, NA, NA),
               col = c("brown3", "brown3", "black", NA),
               fill = c(NA, NA, NA, "grey80"), 
               border = NA )

dev.off()

dat %>%
        mutate(month = month(date)) %>%
        filter(month %in% c(5,6)) %>%
        select(discharge_m3s, temp_C) %>%
        summary()
        # summarize(across(everything(), mean, na.rm = TRUE))
        # summarize(across(everything(), ~mean(., na.rm = TRUE)))

dat <- select(dat, doy, mean_va, p05_va, p95_va, mean_2018) 
dat[,2:5] <- dat[,2:5]*(12*2.54/100)^3  # convert from ft3 to m3

#########################
# Compare baseflow to mean discharge from south Durham WWTP on Farrington Rd. 
# According to the city of durham, the plant in permitted to discharge 20 million gallons per day
# and is currently operating at half capacity, so 10 MGD

wwtp.permit <- 20 #million gallons per day
wwtp.permit <- wwtp.permit *3.78541/1000   # million m3 per day
wwtp.permit <- wwtp.permit*1000000/24/60/60 # m3 per second
wwtp.discharge <- wwtp.permit/2

NHCbaseflow <- quantile(dat$mean_va, 0.2) # baseflow defined as 20th percentile of flow. 
NHCbaseflow2 <- baseflow_sep(dat$mean_va)
NHCbaseflow2018 <- baseflow_sep(dat$mean_2018)#, parms=c(window_size=10, f_lowess=0.1))

plot(dat$mean_2018, type="l", log="y")
lines(NHCbaseflow2018, col="red")
wwtp.flow.contribution <- wwtp.discharge/NHCbaseflow2018[146]*100
wwtp.fc.max <- wwtp.discharge/(min(dat$mean_va, na.rm=T))*100
wwtp.fc.min <- wwtp.discharge/(wwtp.discharge+max(NHCbaseflow2018, na.rm=T))*100
quantile(dat$mean_va, .5, na.rm = T)

 #####################################################

dat.ext <- rbind(dat[364:366,],dat,dat[1:3,])
dat <- data.frame(apply(dat, 2, na.fill, fill="extend"))

dat$p05_va_roll <- rollmean(dat.ext$p05_va, 7)
dat$p95_va_roll <- rollmean(dat.ext$p95_va, 7)

 ggplot(dat, mapping = aes(x=doy, y=p95_va))+
   geom_ribbon(aes(ymin=p05_va, ymax=p95_va))

plot(dat$doy, dat$p95_va, log="y",ylim = c(.1,200), type ="n", 
     ylab = expression(paste("Discharge (m"^3,"/s)",sep="" )), 
     xlab = "day of year", cex.main = 1,
     main = "NHC 2018 hydrograph at Blands USGS station")
polygon(x=c(dat$doy, rev(dat$doy)),y=c((dat$p05_va), rev((dat$p95_va))), 
        col = alpha("grey50",.3), border=NA)

lines(dat$doy, dat$mean_2018, lwd = 1.5, col = "grey30")
d<- format(as.Date("2018-05-25"), '%j')
points(d, dat[dat$doy==d,]$mean_2018, pch = 20, cex=1.5)


baseflow_sep <- function (runoff, method = "DFM", parms = c(c = 0.925, window_size = 10, 
                                            f_lowess = 0.1)) 
{
        lr = length(runoff)
        if (method == "DFM") {
                cc = parms["c"]
                qd = array(runoff[1], lr)
                diff_runoff = diff(runoff)
                qd = c(runoff[1], filter(x = (1 + cc)/2 * diff_runoff, 
                                         filter = cc, method = "recursive", init = runoff[1]))
                if (any(qd < 0)) 
                        for (i in 2:lr) qd[i] = max(0, cc * qd[i - 1] + 
                                                            (1 + cc)/2 * diff_runoff[i - 1], na.rm = TRUE)
                qb = runoff - qd
        }
        if (method == "constant_slope") {
                min_points = which(diff(sign(diff(runoff))) > 0) + 1
                slope = parms["c_slope"]
                if (is.na(slope)) 
                        slope = diff(range(runoff))/lr/2
                qb = runoff
                while (length(min_points) > 0) {
                        i = min_points[1]
                        qb[i:lr] = runoff[i] + (0:(lr - i)) * slope
                        tt = which(runoff < qb)
                        if (any(tt)) {
                                intersect_i = min(tt)
                                qb[intersect_i:lr] = runoff[intersect_i:lr]
                                min_points = min_points[min_points >= intersect_i]
                        }
                        else break
                }
        }
        if (method == "RLSWM") {
                qb = array(0, lr)
                window_size = parms["window_size"]
                f_lowess = parms["f_lowess"]
                for (i in 1:length(runoff)) {
                        recent_q = runoff[max(1, i - window_size):i]
                        qb[i] = min(recent_q)
                }
                base_flow = lowess(qb, f = f_lowess)$y
                ratio = qb/runoff
                if (max(ratio, na.rm=T) > 1) 
                        qb = qb/max(ratio, na.rm=T)
        }
        return(qb)
}
