### Plot average annual hydrograph for USGS Blands Station 
# 2020-02-19

# Updates:
#       2020-03-27 Added section to calculate NHC baseflow and compare to wwtp discharge
#       2020-04-06 Updated with baseflow separation package RHydro
library(zoo)
library(ggplot2)
#install.packages("RHydro", repos="http://R-Forge.R-project.org")
library(RHydro)

dat<- read.csv("data/raw/NHC_Blands_annualdata.csv", header=T)


dat$doy <- seq(1,366)

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

plot(dat$doy, dat$p95_va, log="y",ylim = c(.1,500), type ="n", 
     ylab = expression(paste("Discharge (m"^3,"/s)",sep="" )), 
     xlab = "day of year", cex.main = 1,
     main = "NHC 2018 hydrograph at Blands USGS station")
polygon(x=c(dat$doy, rev(dat$doy)),y=c((dat$p05_va), rev((dat$p95_va))), 
        col = "lightblue", border=NA)

lines(dat$doy, dat$mean_2018, lwd = 2, col = "steelblue")
d<- format(as.Date("2018-05-25"), '%j')
points(d, dat[dat$doy==d,]$mean_2018, pch = 20, cex=2)


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
