#####################
# Calculate DO slopes post storm
# A Carter 2020 Jan 28

library(lubridate)

# Read in DO data and daily site summarys:
DOSiteSummary <- read.csv("data/DailyDOsummary.csv", header=T)
DOdat <- read.csv("data/raw/allSites_20180706.csv", header = T, stringsAsFactors = F)
NHCsites2018 <- c('Mud','MC751','MC3','MC2','MC1','UNHC',
                  'NHC','NHC5','NHC4','NHC3','NHC2','NHC1')
DOdat$DateTime_UTC <- ymd_hms(DOdat$DateTime_UTC)
DOdat$DateTime <- with_tz(DOdat$DateTime_UTC, tz="EST")
# The targeted storm occured on Jun 26:

startdate <- ymd_hms("2018-06-26 13:00:00", tz="EST")
enddate <- ymd_hms("2018-06-29 18:00:00", tz="EST")
plot(DOdat[DOdat$site=="MC3",]$DateTime, DOdat[DOdat$site=="MC3",]$DO_mgL, 
     xlim = c(startdate,enddate))

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

    