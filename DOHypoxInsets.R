###############

library(lubridate)
library(tidyverse)
library(zoo)

# Read in compiled DO data
DOdat <- read_csv("data/raw/allSites_20180706.csv")
DOdat$DateTime <- with_tz(DOdat$DateTime_UTC, tz="EST")
NHCsites2018 <- c("Mud","MC751","MC3","MC2","MC1","UNHC","NHC",
                  "NHC5","NHC4","NHC3","NHC2","NHC1")
start_date <- ymd_hms("2018-06-12 13:45:00 UTC")
end_date <- ymd_hms("2018-07-06 02:45:00 UTC")
DOdat$site <- factor(DOdat$site, levels = NHCsites2018)


# Make plot of DO with an inset bar of hypoxia
DOsummaryStats <- read_csv("data/DOtimeseriesSummary.csv")

ggplot(DOdat, aes(DateTime, DO_mgL, group = 1)) +
  geom_line(color="black")+
  facet_wrap(~site)

png("figures/DOhypoxiaInsets.png", width=6.5, height=4, units="in", res=300)
par(mfrow = c(3,4), cex = 1,
    mar = c(0,0,0,0), oma = c(3,2.5,1.5,2.5))
for(i in 1:length(NHCsites2018)){
  ts = DOdat[DOdat$site==NHCsites2018[i],]$DateTime
  DO = 100*DOdat[DOdat$site==NHCsites2018[i],]$persatDO 
  plot(ts,DO, xaxt='n', yaxt='n', type = 'l', ylim = c(0,120))
   # plot level data in light grey: 
  level <-DOdat[DOdat$site==NHCsites2018[i],]$Level_m
  if(!is.na(level)){
    par(new=T)
    plot(ts,level,type="n", axes=FALSE, ylim = c(0,5))
    polygon(x=c(ts, rev(ts)), y=na.approx(c(level, rep(-1,length(level)))),
            col = "slategray3",  border=NA)
    if(i %in% c(4,8,12)){
      axis(4, at=c(0,2,4),col = "grey30", col.axis = "grey20", tck=-.05, labels=NA)
      axis(4, at=c(0,2,4),col = "grey30",lwd = 0, line = -.6, cex.axis=.7)
      }
    par(new=T)
    plot(ts,DO,xaxt="n",yaxt="n",xlab="", ylab="", type = 'l', ylim = c(0,120))
    
  }
  # plot storm date
  abline(v=ymd_hms("2018-06-26 16:00:00"), col="grey60", lty=2)

  if(NHCsites2018[i]=="MC751"){
    text(x= ymd_hms("2018-06-29 00:00:00"), y=110, labels="MC4", adj = c(.15,0), cex=.9)
  } else if(NHCsites2018[i]=="Mud"){
    text(x= ymd_hms("2018-06-29 00:00:00"), y=110, labels="Mtrib", adj = c(.15,0), cex=.9)
  } else text(x= ymd_hms("2018-06-29 00:00:00"), y=110, labels=NHCsites2018[i], adj = c(.15,0), cex=.9)
  startbox<- ymd_hms("2018-06-13 00:00:00")
  endbox <- ymd_hms("2018-06-22 00:00:00")
  boxlength <- endbox-startbox
  x <- c(startbox,endbox,endbox,startbox)
  y <- c(106, 106, 120,120)
  polygon(x=x, y = y, col = "white")
  stats <- DOsummaryStats[DOsummaryStats$site==NHCsites2018[i],]
  x50 <- stats$prHypox50*boxlength
  x50 <- c(startbox, startbox+x50, startbox+x50, startbox)
  polygon(x=x50, y=y, col = "grey80")
  # x3 <- stats$prHypox3*boxlength
  # x3 <- c(startbox, startbox+x3, startbox+x3, startbox)
  # polygon(x=x3, y=y, col = "grey40")
  x0 <- stats$prHypox0*boxlength
  x0 <- c(startbox, startbox+x0, startbox+x0, startbox)
  polygon(x=x0, y=y, col = "grey20")
  # Add axes
  if(i %in% c(1,5,9)){
    axis(2, at=c(0,50,100),col = "grey30", col.axis = "grey20", tck=-.05, labels=NA)
    axis(2, at=c(0,50,100),col = "grey30",lwd = 0, line = -.6, cex.axis=.7)
  }
  
  if(i %in% 9:12){
    t <- seq(start_date+60*60*24*4, end_date,by = "8 days")
    axis(side = 1, at = t, labels=FALSE) 
    text(x = t, y = par("usr")[3]-3 , labels = format(t, "%m-%d"), 
         pos = 1, xpd=NA, col = "grey30", cex=.7)
    
  }
}
mtext(text="Date",side=1,line=1.5,outer=TRUE)
mtext(text="DO (%sat)",side=2,line=1.5,outer=TRUE)
mtext(text="level (m)",side=4,line=1.5,outer=TRUE)
par(new=T, mfrow = c(1,1), mar = c(0,0,0,0), oma = c(0,0,0,0))

legend("top", legend=c("Storm event", "Water level ", "pr(DO% sat < 50)   ", "pr(DO% sat = 0)    "), 
       cex=.8, bty="n",ncol=4, xpd=NA,
       col=c("grey50",NA,NA,NA),
       lty=c(2,NA,NA,NA),
       fill=c(NA,"slategray3", "grey80", "grey20"),
       border=c(NA,"slategray3","black","black"),
       x.intersp=c(1,.2,.2,.2))

dev.off()
