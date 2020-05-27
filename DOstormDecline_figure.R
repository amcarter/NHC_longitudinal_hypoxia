#####################
# Calculate DO slopes post storm
# A Carter 2020 Jan 28
# Update 2020 May 22
#   cleaned up to be more generic, run on more cases

library(lubridate)
library(dplyr)

setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_longitudinal_hypoxia/Code and Figs")
# Read in DO data and daily site summarys:
DOSiteSummary <- read.csv("data/DailyDOsummary.csv", header=T)
DOdat <- read.csv("data/raw/allSites_20180706.csv", header = T, stringsAsFactors = F)
NHCsites2018 <- c('Mud','MC751','MC3','MC2','MC1','UNHC',
                  'NHC','NHC5','NHC4','NHC3','NHC2','NHC1')
DOdat$DateTime_UTC <- ymd_hms(DOdat$DateTime_UTC)
DOdat$DateTime <- with_tz(DOdat$DateTime_UTC, tz="EST")
# The targeted storm occured on Jun 26:
stormdate <- as.Date("2018-06-26")
startdate <- ymd_hms("2018-06-25 00:00:00 EST")
enddate <- ymd_hms("2018-07-02 00:00:00 EST")


# Set netagive DO values to zero. Do I need to refine this?
w <-which(DOdat$DO_mgL<0)
DOdat[w, c(3,7)] <- 0

DOdat$site <- factor(DOdat$site, levels = NHCsites2018)

png("figures/DOStormTrajectories.png", width=6.5, height=4, units="in", res=300)
par(mfrow=c(3,4), mar=c(0,0,0,0), oma=c(5,5,2,2))
for(i in 1:12){
    minwin=0:6
    maxwin=9:19
    rec_days=5
  if(NHCsites2018[i]=="Mud"){
    minwin=c(17,18,19,20,21,22,23,0)
    maxwin=5:16
  }
  if(NHCsites2018[i] %in% c("MC2","MC3", "NHC", "NHC3")){
    rec_days=4
  }  
  tmp<- DOdat[DOdat$site==NHCsites2018[i],]
  
  a<- calc_DO_storm_recovery(tmp$DateTime, tmp$persatDO*100, stormdate, minwin, maxwin, rec_days)
  mtext(paste0(tmp$site[1], "     dmin/day = ",round(-a$params$dDO.day_min,0), "%   "),
        side=3,line=-1.1, adj=1, cex=.5 )
  mtext(paste0("damp/day = ", round(a$params$amp_recovery_percent.day, 0), "%   "), 
        side=3, line=-1.9, adj=1, cex=.5)
  if(i %in% c(1,5,9)){
    axis(2, at=c(0,50,100),col = "grey30", col.axis = "grey20", tck=-.05, labels=NA)
    axis(2, at=c(0,50,100),col = "grey30",lwd = 0, line = -.6, cex.axis=.7)
  }
  
  if(i %in% 9:12){
    t <- seq(startdate+60*60*24*2, enddate-60*60*24,by = "2 days")
    axis(side = 1, at = t, labels=FALSE) 
    text(x = t, y = par("usr")[3]-3 , labels = format(t, "%m-%d"), 
         pos = 1, xpd=NA, col = "grey30", cex=.7)
  }
  mtext(text="Date",side=1,line=1.5,outer=TRUE)
  mtext(text="DO (%sat)",side=2,line=1.5,outer=TRUE)
}

dev.off()


###########################################################
# fraction night hypoxia

par(mfrow=c(3,4), mar=c(0,0,0,0), oma=c(5,5,2,2))
for(i in 1:12){
  tmp<- DOdat[DOdat$site==NHCsites2018[i],]
  a<- calc_night_hypoxia(tmp$DateTime_UTC, tmp$persatDO, .5, 
                         lat=sites$Lat[sites$site==NHCsites2018[2]],
                         long=sites$Long[sites$site==NHCsites2018[2]])
  mtext(paste0(NHCsites2018[i]," ",round(100*a$per_night_hypox,0),"% night   ",round(100*a$per_day_hypox,0),"% day"),1,-1.5)
}





DOdat<- DOdat[which(DOdat$DateTime<=enddate),]
DOdat.mintimes <- DOdat[which(DOdat$Hour %in% c(0,1,2,3,4,5,6)),]
DOdat.mintimes <- DOdat.mintimes[which(DOdat.mintimes$Date >as.Date("2018-06-23")),]
DOdatdaily<- DOdat %>% group_by(site, Date)%>%
  summarize(mean.persatDO = mean(persatDO, na.rm=T),
            sd.persatDO = sd(persatDO, na.rm=T),
            max.persatDO = max(persatDO, na.rm=T),
            time.max.persatDO = DateTime[which(persatDO==max(persatDO, na.rm=T))[1]])
DOdatdaily.mintimes <- DOdat.mintimes %>% group_by(site,Date)%>%
  summarize(min.persatDO = min(persatDO, na.rm=T), 
            time.min.persatDO = DateTime[which(persatDO==min(persatDO, na.rm=T))[ceiling(length(which(persatDO==min(persatDO, na.rm=T)))/2)]])
DOdatdaily <- full_join(DOdatdaily, DOdatdaily.mintimes, by = c("site","Date"))          




png("figures/DOStormDecline.png", width=6.5, height=4, units="in", res=300)
par(mfrow = c(3,4), cex = 1,
    mar = c(0,0,0,0), oma = c(3,2.5,1.5,2.5))
for(i in 1:12){
  # dAmplitude/dtime
  DO <- DOdatdaily[DOdatdaily$site==NHCsites2018[i]&DOdatdaily$Date %in% c(stormDate, stormDate+1, stormDate+2, 
                                                                           stormDate+3, stormDate+4),]
  DO$min.persatDO[DO$Date==stormDate] <- DO$max.persatDO[DO$Date==stormDate]
  DO$time.min.persatDO[DO$Date==stormDate]<- DO$time.max.persatDO[DO$Date==stormDate]
  m<- lm(100*min.persatDO~time.min.persatDO, data=DO)

  ts <- DOdat$DateTime[DOdat$site==NHCsites2018[i]]
  D <- 100*DOdat$persatDO[DOdat$site==NHCsites2018[i]]
  plot(ts,D, xaxt='n', yaxt='n', type = 'l', ylim = c(0,120)) 
  text(x= ymd_hms("2018-06-27 05:00:00"), y=110, labels=paste0(NHCsites2018[i]," ",-round(m$coef[2]*60*60*24),"% dDO/day"),
       adj = c(.15,0), cex=.6)
  points(DO$time.min.persatDO, 100*DO$min.persatDO, pch = 20,col = "brown3")
  abline(m$coef[1], m$coef[2], col = "brown3", lty=2)


  # Add axes
  if(i %in% c(1,5,9)){
    axis(2, at=c(0,50,100),col = "grey30", col.axis = "grey20", tck=-.05, labels=NA)
    axis(2, at=c(0,50,100),col = "grey30",lwd = 0, line = -.6, cex.axis=.7)
  }
  
  if(i %in% 9:12){
    t <- seq(startdate+60*60*24*2, enddate-60*60*24,by = "2 days")
    axis(side = 1, at = t, labels=FALSE) 
    text(x = t, y = par("usr")[3]-3 , labels = format(t, "%m-%d"), 
         pos = 1, xpd=NA, col = "grey30", cex=.7)
    
  }
}
mtext(text="Date",side=1,line=1.5,outer=TRUE)
mtext(text="DO (%sat)",side=2,line=1.5,outer=TRUE)

dev.off()



DOsummarystats <- function(site, DOdat){
  daily <- DOdat[,c("date",paste0(site,".DO_mgL"))] %>% 
    rename("DO"=paste0(site,".DO_mgL"))%>%
    group_by(date) %>%                  
    summarize(mean =mean(DO,na.rm = T),
              min = min(DO, na.rm = T),
              max = max(DO,na.rm = T),
              amp = (max-min))
}
# Gather DO slope points by selecting the storm peak 
# then the min for each of the three days following



dat <- DOdat[,c("site","DateTime", "DO_mgL")]
dat<- dat[-(which(dat$DateTime<startdate)),]
dat<- dat[-(which(dat$DateTime>enddate)),]

dat$site <- factor(dat$site, levels = NHCsites2018)
timesteps <- c(startdate, seq(ymd_hms("2018-06-26 18:00:00", tz="EST"),
                                enddate, by = "days"))
dat_group <- dat %>% 
  mutate(fac=cut(DateTime, timesteps, labels=letters[1:4]))

dat_points<- data.frame() 

for(i in 1:length(NHCsites2018)){
  site <- NHCsites2018[i]
  datsite <- dat_group[dat_group$site==site,]
  tmp <- datsite[datsite$fac==letters[1],]
  DO <-  max(tmp$DO_mgL, na.rm=T) 
  DateTime <- tmp$DateTime[which.max(tmp$DO_mgL)]
  dat_points <- dplyr::bind_rows(dat_points, data.frame(site=site, DO_mgL=DO, DateTime=DateTime))
  for(j in 2:4){
    tmp <- datsite[datsite$fac==letters[j],]
    DO <- min(tmp$DO_mgL, na.rm=T) 
    DateTime <- tmp$DateTime[which.min(tmp$DO_mgL)]
    dat_points <- dplyr::bind_rows(dat_points, data.frame(site=site, DO_mgL=DO, DateTime=DateTime))
  }
}


plot(DOdat[DOdat$site=="MC3",]$DateTime, DOdat[DOdat$site=="MC3",]$DO_mgL, 
     xlim = c(startdate,enddate))

#Calculate slope of decline for each site, save as EOD:
EOD <- data.frame()
for(i in 1:length(NHCsites2018)){
  site <- NHCsites2018[i]
  model <- lm(dat_points[dat_points$site==site,]$DO_mgL~dat_points[dat_points$site==site,]$DateTime)
  resp <- -model$coef[2]*60*60*24
  EOD <- dplyr::bind_rows(EOD, data.frame(site=site, EOD.mgLday=resp))
}

# Add this column to the DO statistics file:

DOmetrics <- read.csv("data/DOtimeseriesMetrics.csv", header=T)
DOmetrics <- full_join(DOmetrics, EOD, by = "site")

write.csv(DOmetrics, file = "data/DOtimeseriesMetrics.csv", row.names = F)


PlotDOdeclinefit <- function(site,DOdat, declines){
    plot(DOdat[DOdat$site==site,]$DateTime, DOdat[DOdat$site==site,]$DO_mgL,
       xlim = c(startdate,enddate), ylim = c(0,10), 
       type="l", lwd = 1.5, xaxt="n", yaxt="n")
    points(declines[declines$site==site,]$DateTime, declines[declines$site==site,]$DO_mgL,
                    col = "red", pch = 20, cex = 1.2)
    model <- lm(declines[declines$site==site,]$DO_mgL~declines[declines$site==site,]$DateTime)
    par(new=T)
    abline(model$coef[1], model$coef[2], lty=2, col = "red")
    resp <- -model$coef[2]*60*60*24
    text((startdate+20000), 10, adj=c(0,1), labels = site)
    text((startdate+120000), 10, adj=c(0,1), labels = paste0(round(resp, 1), " mg/L/d O2"))
}

par(mfrow = c(3,4), cex = 1,
    mar = c(0,0,0,0), oma = c(3,3,.5,.5))

for(i in 1:length(NHCsites2018)){
  PlotDOdeclinefit(NHCsites2018[i], DOdat, dat_points)
  if(i %in% c(1,5,9)){ axis(2)}
  if(i %in% c(9,10,11,12)){
    axis(1, lab=c("6-27", "6-28","6-29"), 
         at=seq(ymd_hms("2018-06-27 00:00:00"), by="day", length.out=3))
    }
}



PlotDOdeclinefit("MC751", DOdat, dat_points)    
axis(2)  
PlotDOdeclinefit("MC3", MC3)    
PlotDOdeclinefit("MC2", MC2)    
PlotDOdeclinefit("MC1", MC1)
axis(2)
PlotDOdeclinefit("NHC5", NHC5)    
PlotDOdeclinefit("NHC4", NHC4)    
PlotDOdeclinefit("NHC3", NHC3)
axis(2)
PlotDOdeclinefit("NHC2", NHC2)    
PlotDOdeclinefit("NHC1", NHC1)    




dat <- DOdat[,c("date","DateTime_UTC", paste0(NHCsites2018[i],".DO_mgL"))] %>% 
  rename("DO"=paste0(NHCsites2018[i],".DO_mgL"))
  tmp <- dat%>% filter(date==date("2018-06-26"))
  MC751 <- tmp[which.max(tmp$DO),]
  tmp <- dat%>% filter(date==date("2018-06-27"))
  MC751 <- rbind(MC751, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-28"))
  MC751 <- rbind(MC751, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-29"))
  MC751 <- rbind(MC751, tmp[which.min(tmp$DO),])
  
  dat <- DOdat[,c("date","DateTime_UTC", paste0(NHCsites2018[2],".DO_mgL"))] %>% 
  rename("DO"=paste0(NHCsites2018[2],".DO_mgL"))
  tmp <- dat%>% filter(date==date("2018-06-26"))
  MC3 <- tmp[which.max(tmp$DO),]
  tmp <- dat%>% filter(date==date("2018-06-27"))
  MC3 <- rbind(MC3, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-28"))
  MC3 <- rbind(MC3, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-29"))
  MC3 <- rbind(MC3, tmp[which.min(tmp$DO),])
  
  
  dat <- DOdat[,c("date","DateTime_UTC", paste0(NHCsites2018[4],".DO_mgL"))] %>% 
    rename("DO"=paste0(NHCsites2018[4],".DO_mgL"))
  tmp <- dat%>% filter(date==date("2018-06-26"))
  MC1 <- tmp[which.max(tmp$DO),]
  tmp <- dat%>% filter(date==date("2018-06-27"))
  MC1 <- rbind(MC1, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-28"))
  MC1 <- rbind(MC1, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-29"))
  MC1 <- rbind(MC1, tmp[which.min(tmp$DO),])
 
   dat <- DOdat[,c("date","DateTime_UTC", paste0(NHCsites2018[3],".DO_mgL"))] %>% 
    rename("DO"=paste0(NHCsites2018[3],".DO_mgL"))
  tmp <- dat%>% filter(date==date("2018-06-26"))
  MC2 <- tmp[which.max(tmp$DO),]
  tmp <- dat%>% filter(date==date("2018-06-27"))
  MC2 <- rbind(MC2, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-28"))
  MC2 <- rbind(MC2, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-29"))
  MC2 <- rbind(MC2, tmp[which.min(tmp$DO),])
  
  dat <- DOdat[,c("date","DateTime_UTC", paste0(NHCsites2018[5],".DO_mgL"))] %>% 
    rename("DO"=paste0(NHCsites2018[5],".DO_mgL"))
  tmp <- dat%>% filter(date==date("2018-06-26"))
  NHC5 <- tmp[which.max(tmp$DO),]
  tmp <- dat%>% filter(date==date("2018-06-27"))
  NHC5 <- rbind(NHC5, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-28"))
  NHC5 <- rbind(NHC5, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-29"))
  NHC5 <- rbind(NHC5, tmp[which.min(tmp$DO),])
  
  
  dat <- DOdat[,c("date","DateTime_UTC", paste0(NHCsites2018[6],".DO_mgL"))] %>% 
    rename("DO"=paste0(NHCsites2018[6],".DO_mgL"))
  tmp <- dat%>% filter(date==date("2018-06-26"))
  NHC4 <- tmp[which.max(tmp$DO),]
  tmp <- dat%>% filter(date==date("2018-06-27"))
  NHC4 <- rbind(NHC4, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-28"))
  NHC4 <- rbind(NHC4, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-29"))
  NHC4 <- rbind(NHC4, tmp[which.min(tmp$DO),])

  
  dat <- DOdat[,c("date","DateTime_UTC", paste0(NHCsites2018[7],".DO_mgL"))] %>% 
    rename("DO"=paste0(NHCsites2018[7],".DO_mgL"))
  tmp <- dat%>% filter(date==date("2018-06-26"))
  NHC3 <- tmp[which.max(tmp$DO),]
  tmp <- dat%>% filter(date==date("2018-06-27"))
  NHC3 <- rbind(NHC3, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-28"))
  NHC3 <- rbind(NHC3, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-29"))
  NHC3 <- rbind(NHC3, tmp[which.min(tmp$DO),])
  
  dat <- DOdat[,c("date","DateTime_UTC", paste0(NHCsites2018[8],".DO_mgL"))] %>% 
    rename("DO"=paste0(NHCsites2018[8],".DO_mgL"))
  tmp <- dat%>% filter(date==date("2018-06-26"))
  NHC2 <- tmp[which.max(tmp$DO),]
  tmp <- dat%>% filter(date==date("2018-06-27"))
  NHC2 <- rbind(NHC2, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-28"))
  NHC2 <- rbind(NHC2, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-29"))
  NHC2 <- rbind(NHC2, tmp[which.min(tmp$DO),])
  
  dat <- DOdat[,c("date","DateTime_UTC", paste0(NHCsites2018[9],".DO_mgL"))] %>% 
    rename("DO"=paste0(NHCsites2018[9],".DO_mgL"))
  tmp <- dat%>% filter(date==date("2018-06-26"))
  NHC1 <- tmp[which.max(tmp$DO),]
  tmp <- dat%>% filter(date==date("2018-06-27"))
  NHC1 <- rbind(NHC1, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-28"))
  NHC1 <- rbind(NHC1, tmp[which.min(tmp$DO),])
  tmp <- dat%>% filter(date==date("2018-06-29"))
  NHC1 <- rbind(NHC1, tmp[which.min(tmp$DO),])
  
  
  sdat <- NHC1

    