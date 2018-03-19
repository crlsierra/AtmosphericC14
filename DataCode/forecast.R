################################################################################
## Code to reproduce all results from
# Forecasting atmospheric radiocarbon decline to pre-bomb values  
################################################################################
#################################################### by Carlos A. Sierra       

#Load requiered R packages
library(SoilR) # The Hua et al. (2013) series are stored in this package
library(forecast)
library(xtable) # Use only to produce latex code from data.frames and create tables

#Load other auxiliary data provided with the manuscript
cities=read.csv("~/SOIL-R/Manuscripts/AtmosphericC14/DataCode/S1_Cities.csv")
Levin=read.table("~/SOIL-R/Manuscripts/AtmosphericC14/DataCode/S3_Levin.txt",header=TRUE)

lambda=1/8267

# Function to transform dataseries to F' from F
FpfromF=function(X){ # X must be a data.frame with first column time and second column Fraction Modern
  y=X[,1]
  fM=X[,2]
  Fp=fM*exp((1950-y)*lambda)
  return(data.frame(year=y,Fp=Fp))
}


Fpnh=FpfromF(Hua2013$NHZone1[,c(1,4)]) #Absolute fraction modern F'
Fpsh=FpfromF(Hua2013$SHZone3[,c(1,4)]) #Absolute fraction modern F'

# Transform the available time series to regularly spaced intervals
y0=1975
m=12 #Monthly time-scale
yNm=seq(y0,2009.6,by=1/m)
nz1=spline(Fpnh,xout=yNm)
ySm=seq(y0,2011.2,by=1/m)
sz3=spline(Fpsh,xout=ySm)

nh=ts(nz1$y,start=y0,freq=m)
sh=ts(sz3$y,start=y0,freq=m)


#### Figure 1
#pdf("~/SOIL-R/Manuscripts/AtmosphericC14/Figures/HuaSeries.pdf")
par(mfrow=c(2,1),mar=c(4,4,1,4))
plot(FpfromF(Hua2013$NHZone1[,c(1,4)]),
     type="l",xlab="Year",ylab="Absolute fraction modern",bty="n")
lines(FpfromF(Hua2013$NHZone2[,c(1,4)]),col=2)
lines(FpfromF(Hua2013$NHZone3[,c(1,4)]),col=3)
lines(FpfromF(Hua2013$SHZone12[,c(1,4)]),col=4)
lines(FpfromF(Hua2013$SHZone3[,c(1,4)]),col=5)
abline(h=1,lty=2)
axis(4,at=seq(1.0,2.0,by=0.2),labels=as.character((seq(1.0,2.0,by=0.2)-1)*1000))
mtext(expression(paste(Delta^14,"C ","(per mil)")),side=4,line=3)
legend("topright",c(
  "Northern hemisphere zone 1",
  "Northern hemisphere zone 2",
  "Northern hemisphere zone 3",
  "Southern hemisphere zones 1 and 2",
  "Southern Hemispher zone 3"
),lty=1,col=1:5,bty="n")
legend("topleft","a",bty="n")

plot(nh,bty="n",ylim=c(1,1.5),ylab="Absolute fraction modern",type="l",xlab="Year AD",col="brown")
lines(sh,col="darkgreen")
abline(h=1,lty=2)
legend("topleft","b",bty="n")
legend("topright",c("Northern hemisphere", "Southern hemisphere"),col=c("brown","darkgreen"),lty=1,bty="n")
axis(4,at=seq(1.0,1.5,by=0.2),labels=as.character((seq(1.0,1.5,by=0.2)-1)*1000))
mtext(expression(paste(Delta^14,"C ","(per mil)")),side=4,line=3)
par(mfrow=c(1,1))
#dev.off()

##################################
# Time series decomposition at the monthly time-scale
mnhm=ets(nh) #Best model selected by AIC
mnhm$method
summary(mnhm)

mshm=ets(sh) #Best model selected by AIC
mshm$method
summary(mshm)

# Observed vs predicted. Not in manuscript
plot(nh,mnhm$fitted, pch=20, xlab="Observations (Absolute fraction modern)", ylab="Predictions (Absolute fraction modern)", ylim=c(1,1.5), xlim=c(1,1.5), col="brown")
points(sh, mshm$fitted, pch=20, col="darkgreen")
abline(a=0,b=1)
legend("topleft", c("Northern hemisphere", "Southern hemisphere"), pch=20, col=c("brown", "darkgreen"), bty="n")

#### Slopes
bnhm=mnhm$states[,"b"]
bshm=mshm$states[,"b"]

#annual slopes
Snh=(tapply(bnhm, floor(time(bnhm)), FUN=sum))*1000
Ssh=(tapply(bshm, floor(time(bshm)), FUN=sum))*1000

#annual error
Enh=(tapply(mnhm$residuals, floor(time(mnhm$residuals)), FUN=function(x){sqrt(sum(x^2))}))*1000
Esh=(tapply(mshm$residuals, floor(time(mshm$residuals)), FUN=function(x){sqrt(sum(x^2))}))*1000

## Seasonality
snhm=mnhm$states[,"s1"]
sshm=mshm$states[,"s1"]

# Table 1
lasty=matrix(NA,7,4)
lasty[1:5,1]=Snh[32:36]
lasty[1:5,2]=Enh[31:35]
lasty[,3]=Ssh[32:38]
lasty[,4]=Esh[31:37]
lasty
xtable(as.data.frame(lasty,row.names = as.character(2005:2011))) #uncomment for table in LaTeX format

### Figure 2
#pdf("~/SOIL-R/Manuscripts/AtmosphericC14/Figures/SlopeSeason.pdf")
par(mfrow=c(3,1),mar=c(4,4,1,4))
plot(mnhm$states[,1], ylab="Level in Absolute fraction modern", bty="n", col="brown", xlab="", ylim=c(1,1.4))
lines(mshm$states[,1], col="darkgreen")
axis(4,at=seq(1.0,1.4,by=0.1),labels=as.character((seq(1.0,1.4,by=0.1)-1)*1000))
mtext(expression(paste(Delta^14,"C ","(per mil)")),side=4,line=3, cex=0.7)
legend("topright", c("Northern hemisphere", "Southern hemisphere"), lty=1,col=c("brown", "darkgreen"), bty="n")

plot(bnhm, ylim=c(-0.004,0),ylab="Growth in Absolute fraction modern", xlab="",bty="n", col="brown")
lines(bshm,col="darkgreen")
axis(4,at=seq(-0.004,0,by=0.001),labels=as.character(seq(-4,0,by=1)))
mtext(expression(paste(Delta^14,"C "," (per mil ",month^-1,")")),side=4,line=3,cex=0.7)

plot(snhm,ylab="Seasonality (unitless)",xlab="Year",bty="n", col="brown",ylim=c(0.996,1.004))
lines(sshm, col="darkgreen")
abline(h=1,lty=3)
par(mfrow=c(1,1))
#dev.off()

# Full decomposition
plot(cbind(observed=nh,level=mnhm$states[,1], slope=mnhm$states[,"b"], season=mnhm$states[,"s1"]),
     main="Northern hemisphere")

plot(cbind(observed=sh,level=mshm$states[,1], slope=mshm$states[,"b"], season=mshm$states[,"s1"]),
     main="Southern hemisphere")

#Level + slope
tmnhm=mnhm$states[,1]+bnhm
tmshm=mshm$states[,1]+bshm
plot(tmnhm, ylab="Trend (Level plus slope)", bty="n", col="brown")
lines(tmshm, col="darkgreen")

# Residuals
#pdf("~/SOIL-R/Manuscripts/AtmosphericC14/Figures/Error.pdf")
par(mfrow=c(2,2), mar=c(4,4,1,1))
plot(mnhm$residuals, ylab="Residuals", col="brown", ylim = c(-0.02, 0.02), xlim = c(1975, 2010), bty="n")
abline(h=0, lty=2)

hist(mnhm$residuals, xlab="Residuals", main="", xlim = c(-0.02, 0.02))

plot(mshm$residuals, ylab="Residuals", col="darkgreen", ylim = c(-0.02, 0.02), xlim = c(1975, 2010), bty="n")
abline(h=0, lty=2)

hist(mshm$residuals, xlab="Residuals", main="", xlim = c(-0.02, 0.02))
par(mfrow=c(1,1))
#dev.off()

###################################
# Forecast
pq=4 #quaterly time scale
yNq=seq(y0,2009.6,by=1/pq)
nz1q=spline(Fpnh,xout=yNq)
ySq=seq(y0,2011.2,by=1/pq)
sz3q=spline(Fpsh,xout=ySq)

nhq=ts(nz1q$y,start=y0,freq=pq)
shq=ts(sz3q$y,start=y0,freq=pq)

# Fit ETS model to quaterly series
mnhq=ets(nhq)
mshq=ets(shq)

foryrs=20

#Use the ETS model to forecast
fnh=forecast(mnhq,h=foryrs*pq, level=c(68, 95))
fsh=forecast(mshq,h=foryrs*pq, level=c(68, 95))

#Fit all Exponential Smoothing Models to European time-series
errortype <- c("A","M")
trendtype <- c("N","A","M")
seasontype <- c("N","A","M")

models=expand.grid(errortype,trendtype,seasontype) #all combinations of possible models. Output is a matrix
models=models[-c(5,11,12,13,15,17),]
modeltype=paste(models[,1],models[,2],models[,3],sep="") #model types concatened as single char string.

# Take independent observations by Levin and calculates the MSR from the forecast.
for(i in 1:length(modeltype)){
  candidate=ets(nhq,model=modeltype[i],damped=NULL)
  forcCand=forecast(candidate,h=foryrs*pq, level=c(68, 95))
  meanForc=forcCand$mean
  forcFunc=splinefun(x=as.numeric(time(meanForc)),y=as.numeric(meanForc))

  obs=(Levin[,5]/1000)+1
  pred=forcFunc(Levin[,1])
  msr=(sum((pred-obs)^2))/length(obs)
  models[i,4]=msr
}

#Selects the model type that best fits the observations for NH.
mEU=ets(nhq,model=modeltype[which.min(models[,4])])
fEU=forecast(mEU, h=foryrs*pq, level=c(69,95))


### Figure 3
#pdf("~/SOIL-R/Manuscripts/AtmosphericC14/Figures/Forecast.pdf")
par(mfrow=c(2,1),mar=c(4,4,1,4))
plot(fnh, main="Northern hemisphere",xlim=c(2000,2030),ylim=c(0.9,1.2),bty="n", ylab="Absolute fraction modern")
axis(4,at=seq(0.9,1.2,by=0.05),labels=as.character((seq(0.9,1.2,by=0.05)-1)*1000))
mtext(expression(paste(Delta^14,"C ","(per mil)")),side=4,line=3,cex=0.7)
abline(h=1)
# points(Levin$Year,(Levin$Schauinsland/1000)+1,pch=20,cex=0.6)
# points(Levin$Year,(Levin$Jungfraujoch/1000)+1,pch=20,col=2,cex=0.6)
#legend("bottomleft",legend=c("Jungfraujoch","Schauinsland"),bty="n",pch=20,col=c(2,1))

plot(fsh, main="Southern hemisphere",xlim=c(2000,2030),ylim=c(0.9,1.2), xlab="Year",bty="n", ylab="Absolute fraction modern")
axis(4,at=seq(0.9,1.2,by=0.05),labels=as.character((seq(0.9,1.2,by=0.05)-1)*1000))
mtext(expression(paste(Delta^14,"C ","(per mil)")),side=4,line=3,cex=0.7)
abline(h=1)
par(mfrow=c(1,1))
#dev.off()

#pdf("~/SOIL-R/Manuscripts/AtmosphericC14/Figures/ForecastEurope.pdf")
par(mfrow=c(2,1),mar=c(4,4,1,4))
plot(fnh, main="Northern hemisphere",xlim=c(2010,2020),ylim=c(0.9,1.2),bty="n", ylab="Absolute fraction modern")
axis(4,at=seq(0.9,1.2,by=0.05),labels=as.character((seq(0.9,1.2,by=0.05)-1)*1000))
mtext(expression(paste(Delta^14,"C ","(per mil)")),side=4,line=3,cex=0.7)
abline(h=1)
points(Levin$Year,(Levin$Schauinsland/1000)+1,pch=20,cex=0.6)
points(Levin$Year,(Levin$Jungfraujoch/1000)+1,pch=20,col=3,cex=0.6)
legend("bottomleft",legend=c("Jungfraujoch","Schauinsland"),bty="n",pch=20,col=c(3,1))
legend("topleft", "a", bty="n")

plot(fEU, main="Central Europe",xlim=c(2010,2020),ylim=c(0.9,1.2),bty="n", ylab="Absolute fraction modern")
axis(4,at=seq(0.9,1.2,by=0.05),labels=as.character((seq(0.9,1.2,by=0.05)-1)*1000))
mtext(expression(paste(Delta^14,"C ","(per mil)")),side=4,line=3,cex=0.7)
abline(h=1)
points(Levin$Year,(Levin$Schauinsland/1000)+1,pch=20,cex=0.6)
points(Levin$Year,(Levin$Jungfraujoch/1000)+1,pch=20,col=3,cex=0.6)
legend("topleft", "b", bty="n")
par(mfrow=c(1,1))
#dev.off()

########
# Calculate the probability of D14 <= 0
probfor=function(m,h,pseq){
  f=forecast(m,h=h,level=pseq)
  p=which(f$lower<=1,arr.ind=TRUE)
  tseq=seq(tsp(f$mean)[1],tsp(f$mean)[2],by=1/tsp(f$mean)[3])
  P=cbind(p,year=tseq[p[,1]],P=100-pseq[p[,2]])
  cumP=aggregate(P[,4],by=list(P[,3]),max)
  return(cumP)
}

foryrs=30

Pnh=probfor(m=mnhq,h=foryrs*pq,pseq=seq(1,99,by=1))
Psh=probfor(m=mshq,h=foryrs*pq,pseq=seq(1,99,by=1))
PEU=probfor(m=mEU,h=foryrs*pq,pseq=seq(1,99,by=1))

# Figure 5
#pdf("~/SOIL-R/Manuscripts/AtmosphericC14/Figures/Prob.pdf")
plot(Pnh,type="l",xlim=c(2010,2040),ylim=c(0,80),xlab="Year",ylab="Probability (%)",bty="n", col="brown")
lines(Psh,col="darkgreen")
lines(PEU,col=4)
abline(h=20,lty=2)
abline(h=5,lty=3)
abline(v=2018.5, lty=3)
legend("topright", c("Northern hemisphere", "Southern hemisphere", "Central Europe"), lty=1, col=c("brown", "darkgreen", 4), bty="n")
#dev.off()

######### 
# Plot forecast with available data

#Figure 6
colpal=heat.colors(16)
#pdf("~/SOIL-R/Manuscripts/AtmosphericC14/Figures/Cities.pdf")
par(mar=c(5,4,1,5))
plot(fnh,xlab="Year", ylab="Absolute fraction modern", ylim=c(0.8,1.6), xlim=c(2012,2018),main=" ",bty="n")
abline(h=1)
points(cities$Year+cities$DOY/365,1+(cities$D14C)/1000,pch=20,col=colpal)
# points(Levin$Year,(Levin$Schauinsland/1000)+1,pch=20,cex=0.6)
# points(Levin$Year,(Levin$Jungfraujoch/1000)+1,pch=20,col=3,cex=0.6)
legend("toprigh",legend=as.character(cities$City),bty="n",pch=20,col=colpal)
# legend("topleft",legend=c("Jungfraujoch","Schauinsland"),bty="n",pch=20,col=c(3,1))
axis(4,at=seq(0.8,1.6,by=0.2),labels=as.character((seq(0.8,1.6,by=0.2)-1)*1000))
mtext(expression(paste(Delta^14,"C ","(per mil)")),side=4,line=3)
#dev.off()

# Table 2
cts=data.frame(LabID=cities$LabID, 
               City=cities$City, Country=c("Colombia", "USA", "Sweden", "Colombia", "USA", "Austria", "Israel", "Germany", "Germany", "Czech Republic", "Sweden", "Sweden", "Sweden", "Colombia", "USA", "Austria"),
               Sampling_Date= paste(cities$Year,"-",cities$Month, "-", cities$Day, sep=""),Delta14C=cities$D14C, sd=cities$D14Csd, Latitude=cities$Latitude, Longitude=cities$Longitude)

cts
#xtable(cts)

# Save forecasts as comma separated files
write.csv(fnh, "~/SOIL-R/Manuscripts/AtmosphericC14/DataCode/forecastNH.csv")
write.csv(fsh, "~/SOIL-R/Manuscripts/AtmosphericC14/DataCode/forecastSH.csv")
write.csv(fEU, "~/SOIL-R/Manuscripts/AtmosphericC14/DataCode/forecastEU.csv")

################################################################################
# End. Last run 12-2017 on R version 3.4.3


