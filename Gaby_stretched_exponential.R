# -------------------------------------------------------------------
# Stretched_Exp_DOC.m: Program to fit stretched exponentials to DOC 
#                      time series.
# 
# Date: 2020 06 15 
# Gabriel Katul and Alice Carter
#  
# Version: 1
# Translated from matlab to R by A carter
#
# Reference: Laherrere, J., and D. Sornette,1998, 
# Stretched exponential distributions in nature and economy,
#“fat tails” with characteristic scales. 
#The European Physical Journal B-Condensed Matter and Complex 
#Systems 2.4: 525-539.
#
#-------------------------------------------------------------------
library(readr)
library(fitdistrplus)
library(lubridate)

setwd(hypox_projdir)

#######################################
#read in data and metadata
sites <- read_csv("NHC_map/NC_synopticSamplingSites.csv")
annual_dat <- read_csv("data/raw/2019SPsites.csv")
annual_dat$DateTime <- with_tz(annual_dat$DateTime_UTC, tz="EST")
annual_dat$date <- as.Date(annual_dat$DateTime, tz = "EST")
annual_dat$year <- year(annual_dat$date)
annual_dat<- annual_dat[annual_dat$year %in% 2017:2019,]
annual_dat$site_year <- paste(annual_dat$site,year(annual_dat$date), sep="_")

dat <- annual_dat[annual_dat$site_year=="MC751_2019",]
#---- Data format: Date, DOC_Sat (#), DO (mg/L), Water T (C), Site-year 
C_sat<- dat$persatDO
DOC <- dat$DO_mgL
Tw <- dat$WaterTemp_C

#----- Create a vector of time
M<- length(DOC)
t=c(0:M-1)*(15/60)*(1/24)

#----- Compare and plot pdfs


DOC_n=DOC[DOC>2*10^-16]
DOC_n <- DOC_n[!is.na(DOC_n)]
yy=hist(DOC_n,100, plot=F)$counts
xx=hist(DOC_n,100, plot=F)$mids
dx=xx[2]-xx[1]
yy=yy/(dx*sum(yy)) # convert to density
wb <- fitdist(DOC_n, "weibull")
parmHat <- wb$estimate
parmCI <- wb$sd
y <- dweibull(xx,parmHat[1],parmHat[2])

par(mfrow=c(2,2))
plot(xx,yy, xlab="DO (mg/L)", ylab = "pdf")
lines(xx,y)

plot(xx,yy,log="xy", xlab="DO (mg/L)", ylab = "pdf")
lines(xx,y)
yyc=dx*cumsum(yy)
yc=dx*cumsum(y)
plot(xx,yyc,xlab="DO (mg/L)", ylab="cdf")
lines(xx,yc)

eyc=1-yyc
ec=1-yc
plot(xx,eyc)
lines (xx,ec)



#----- build pdf for % sat data
pdf(file="weibull_fits.pdf", onefile=T)
for(i in unique(annual_dat$site_year)){
  dat<- annual_dat[annual_dat$site_year==i,]
  C_sat_n <- dat$persatDO
  C_sat_n <- C_sat_n[!is.na(C_sat_n)]
  C_sat<- C_sat_n[C_sat_n>2*10^-16]# <- 2*10^-16
  yy_n <- hist(C_sat_n, 100, plot=F)$density
  xx_n <- hist(C_sat_n, 100, plot=F)$mids
  
  yy=hist(C_sat,100, plot=F)$counts
  xx=hist(C_sat,100, plot=F)$mids
  dx=xx[2]-xx[1]
  yy=yy/(dx*sum(yy)) # convert to density
  wb <- fitdist(C_sat, "weibull")
  parmHat <- wb$estimate
  parmCI <- wb$sd
  y <- dweibull(xx,parmHat[1],parmHat[2])
  
  par(mfrow = c(2,2), mar = c(4,4,1,1), oma=c(1,1,3,1))
  plot(xx,yy, xlab="DO pctsat", ylab = "pdf")
  lines(xx,y)
  
  plot(xx,yy, xlab="DO pctsat", ylab = "pdf", log="xy")
  lines(xx,y)
  
  yyc=dx*cumsum(yy)
  yc=dx*cumsum(y)
  plot(xx,yyc,xlab="DO pctsat", ylab="cdf")
  lines(xx,yc)
  
  eyc=1-yyc
  ec=1-yc
  plot(xx,eyc,xlab="DO pctsat", ylab="exc")
  lines (xx,ec)
  mtext(i, side=3, outer=T, line=.5, cex=1.2)
}  
dev.off()
