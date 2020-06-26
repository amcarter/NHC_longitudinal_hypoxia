#####################
# Conceptual figure for NHC hypoxia paper

####################################################################
# Panel 1 daily curves based on sin function
par(mar = c(3,3,1,1))
plot(1,type = "n", xlim = c(-pi/2,3*pi/2), ylim = c(0,120))
curve(60*sin(x)+60, add = T, lwd = 2)
plot(x,productive, type = "l")
####################################################################
# Panel 2: Storm recovery across sites

startdate <- ymd_hms("2018-06-26 0:00:00", tz="EST")
enddate <- ymd_hms("2018-07-02 24:00:00", tz="EST")
DOdat$persatDO.smooth<- rollmean(DOdat$persatDO, 36, na.pad=T)
# Plot up:
plot(DOdat[DOdat$site==up,]$DateTime, 
     100*DOdat[DOdat$site==up,]$persatDO.smooth, 
     xlim = c(startdate,enddate), ylim = c(0, 120), 
     col = "darkblue", type = "l", lwd = 2,
     ylab = "DO (% sat)", xlab = "Date")
#plot high
for(i in 1:length(high)){
  lines(DOdat[DOdat$site==high[i],]$DateTime, 
        100*DOdat[DOdat$site==high[i],]$persatDO.smooth,
        col = high.c[i], lwd = 2)
}
# Plot low:
for(i in 1:2){
  lines(DOdat[DOdat$site==low[i],]$DateTime, 
        100*DOdat[DOdat$site==low[i],]$persatDO.smooth,
        lty = low.t[i], lwd = 2)
}
# Plot Steep:
for(i in 1:2){
  lines(DOdat[DOdat$site==steep[i],]$DateTime, 
        100*DOdat[DOdat$site==steep[i],]$persatDO.smooth,
        col = "brown3", lty = low.t[i], lwd = 2)
}
lines(DOdat[DOdat$site==up,]$DateTime, 
      100*DOdat[DOdat$site==up,]$persatDO.smooth, 
      xlim = c(startdate,enddate), ylim = c(0, 120), 
      col = "darkblue", type = "l", lwd = 2)

