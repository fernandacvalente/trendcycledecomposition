library(spdep)
library(sp) 
library(raster)  
library(rasterVis) 
library(maptools)
library(rgeos)
library(dismo)
library(rgdal)
library(INLA)
library(lattice)
library(gridExtra)
library(splancs)
library(xtable)
library(zoo)
library(ggplot2)

##########################################################################################################
## Load data
## 

## Border of São Paulo state 
border=shapesp@polygons
bordersp=border[[1]]@Polygons[[1]]@coords

## Fire occurrence in São Paulo state (point process)
require(splancs)
xy.in <- inout(coordinates(firesco), bordersp)
fires=firesco[xy.in,]

n<- length(fires)
locs=coordinates(fires)
bounds <- gUnaryUnion(shape.sp)

##########################################################################################################
## Quarterly aggregation for daily data
month=as.matrix(substr(firesco$ACQ_DATE,1,7))
unique.month=unique(month)

quarter=month
for (j in 1:length(month)){
  if (substr(month[j],6,7)=="01" | substr(month[j],6,7)=="02" |substr(month[j],6,7)=="03") quarter[j]=1 
  if (substr(month[j],6,7)=="04" | substr(month[j],6,7)=="05" |substr(month[j],6,7)=="06") quarter[j]=2 
  if (substr(month[j],6,7)=="07" | substr(month[j],6,7)=="08" |substr(month[j],6,7)=="09") quarter[j]=3 
  if (substr(month[j],6,7)=="10" | substr(month[j],6,7)=="11" |substr(month[j],6,7)=="12") quarter[j]=4 
}

##########################################################################################################
## Time index
time=double(dim(quarter)[1])+1
time[1]=1
t=1
for (j in 2:dim(quarter)[1]){
  if (quarter[j]!=quarter[j-1]) {
    t=t+1
    time[j]=t}
  if (quarter[j]==quarter[j-1]) {
    time[j]=t}
}

##########################################################################################################
## INLA spatial mesh
smesh <- inla.mesh.2d(boundary=inla.sp2segment(bounds),
                      max.edge=1, cutoff=0.3)

##########################################################################################################
## INLA temporal mesh
n.quarter= max(time)
k <- n.quarter
tmesh <- inla.mesh.1d(seq(0, n.quarter, length=k))

##########################################################################################################
## In order to obtain the space-time aggregation, we find to which polygon belongs each data point in the 
## spatial mesh and to which part of the time belongs each data point. Hereafter, we use these both 
## identification index sets to aggregate the data. 
library(deldir)
dd <- deldir(smesh$loc[,1], smesh$loc[,2])
tiles <- tile.list(dd)

polys <- SpatialPolygons(lapply(1:length(tiles), function(i)
{ p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
n <- nrow(p)
Polygons(list(Polygon(p[c(1:n, 1),])), i)
}))

## Number of fire occurrences in each polygon in the spatial mesh
area <- factor(over(SpatialPoints(cbind(locs[,1], locs[,2])),
                    polys), levels=1:length(polys))
otime=time
t.breaks <- sort(c(tmesh$loc[c(1,k)],tmesh$loc[2:k-1]/2 + tmesh$loc[2:k]/2))
table(time <- factor(findInterval(time, t.breaks),levels=1:(length(t.breaks)-1)))
## Number of fire occurrences in each part of the time in the temporal mesh
time <- factor(otime,levels=1:(length(t.breaks)-1))


agg.dat <- as.data.frame(table(area, time))

## Aggregated data
for(j in 1:2) 
  agg.dat[[j]] <- as.integer(as.character(agg.dat[[j]]))
str(agg.dat)

library(rgeos)
sum(w.areas <- sapply(1:length(tiles), function(i)
{ p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
n <- nrow(p)
pl <- SpatialPolygons(list(Polygons(list(Polygon(p[c(1:n, 1),])), i)))
if (gIntersects(pl, bounds))
  return(gArea(gIntersection(pl, bounds)))
else return(0)
}))

summary(w.areas)

gArea(bounds)
(w.t <- diag(inla.mesh.fem(tmesh)$c0))
(i0 <- n / (gArea(bounds) * diff(range(tmesh$loc))))
## The space-time volumn (area unit per time unit) at each polygon and time knot
summary(e0 <- w.areas[agg.dat$area] * (w.t[agg.dat$time]))

##########################################################################################################
## Setting trend, seasonal and cyclic components
rwtime=agg.dat$time
seastime=agg.dat$time
artime=agg.dat$time

##########################################################################################################
## Setting covariates
loci1=smesh$loc[agg.dat$area,1:2]
vtime=agg.dat$time

## Modis land cover classification is an yearly data. Turn it in quarterly data.
index=vtime==1 | vtime==2 | vtime==3 | vtime==4
ct1=loci1[index,]
vsolo2003=extract(tsol2003,ct1)

index=vtime==5 | vtime==6 | vtime==7 | vtime==8
ct1=loci1[index,]
vsolo2004=extract(tsol2004,ct1)

index=vtime==9 | vtime==10 | vtime==11 | vtime==12
ct1=loci1[index,]
vsolo2005=extract(tsol2005,ct1)

index=vtime==13 | vtime==14 | vtime==15 | vtime==16
ct1=loci1[index,]
vsolo2006=extract(tsol2006,ct1)

index=vtime==16 | vtime==18 | vtime==19 | vtime==20
ct1=loci1[index,]
vsolo2007=extract(tsol2007,ct1)

index=vtime==21 | vtime==22 | vtime==23 | vtime==24
ct1=loci1[index,]
vsolo2008=extract(tsol2008,ct1)

index=vtime==25 | vtime==26 | vtime==27 | vtime==28
ct1=loci1[index,]
vsolo2009=extract(tsol2009,ct1)

index=vtime==29 | vtime==30 | vtime==31 | vtime==32
ct1=loci1[index,]
vsolo2010=extract(tsol2010,ct1)

index=vtime==33 | vtime==34 | vtime==35 | vtime==36
ct1=loci1[index,]
vsolo2011=extract(tsol2011,ct1)

index=vtime==37 | vtime==38 | vtime==39 | vtime==40
ct1=loci1[index,]
vsolo2012=extract(tsol2012,ct1)

index=vtime==41 | vtime==42 | vtime==43 | vtime==44
ct1=loci1[index,]
vsolo2013=extract(tsol2013,ct1)

index=vtime==45 | vtime==46 | vtime==47 | vtime==48
ct1=loci1[index,]
vsolo2014=extract(tsol2014,ct1)

index=vtime==49 | vtime==50 | vtime==51 | vtime==52
ct1=loci1[index,]
vsolo2015=extract(tsol2015,ct1)

index=vtime==53 | vtime==54 | vtime==55 | vtime==56
ct1=loci1[index,]
vsolo2016=extract(tsol2016,ct1)

vsolo=c(vsolo2003,vsolo2004,vsolo2005,vsolo2006,vsolo2007,vsolo2008,vsolo2009,vsolo2010,vsolo2011,vsolo2012,vsolo2013,
        vsolo2014,vsolo2015,vsolo2016)
vsolo2=vsolo==12*1

## Köppen climate classification for São Paulo state
flk=raster("Koppen_Brazil_2013/Koppen Brazil 2013/Raster/koppen_paper/w001000.adf")
koppen=extract(flk,loci1)

## Maximum temperature for São Paulo state, using Laurini* (2019) method.
vtemp=c()
for (jj in 1:56){ #56 quarters
  tempbr= raster("rastersudestetempmax.tif",band=168+jj) 
  indextemp=vtime==jj
  ct1=loci1[indextemp,]
  vtempi=extract(tempbr,ct1)
  vtemp=c(vtemp,vtempi)
}

## Rainfall for São Paulo state, using Laurini* (2019) method.
vpluv=c()
for (jj in 1:56){ #56 quarters
  tempbr= raster("rastersudestepluv.tif",band=168+jj) 
  indextemp=vtime==jj
  ct1=loci1[indextemp,]
  vpluvi=extract(tempbr,ct1)
  vpluv=c(vpluv,vpluvi)
}

## Soil slope for São Paulo state
decliv <- raster("declividade_br.asc")
declivsp<-extract(decliv,loci1)
dummydecliv = declivsp 
dummydecliv[]<- ifelse(declivsp[]>12,1,0) #set as a dummy 

## Sugarcane Agroecological Zoning for the production of ethanol and sugar
zoning<- readOGR("ZoneamentoCana_4C.shp")
zone<-extract(zoning,loci1)

## Monthly global future price of sugar (U.S. centsper pound)
library(xts)
mprice<-read.table("sugar11.csv", sep=";", header=T)
mpricet=ts(mprice[,2],frequency = 12, start = c(1990, 1))
mpricetx=as.xts(mpricet)
mpriceq=apply.quarterly(mpricetx, mean)
ps=mpriceq[53:108,]

price = as.matrix(agg.dat$time)
vprice=c()
for (jj in 1:56){
  indexx = which(as.matrix(agg.dat$time)==jj)
  price[indexx,]=as.numeric(ps)[jj]
}

##########################################################################################################
## SPDE model with Matérn covariance structure
spde <- inla.spde2.matern(smesh)

##########################################################################################################
## Projector matrix
A.st <- inla.spde.make.A(smesh, smesh$loc[agg.dat$area,],group=agg.dat$time, mesh.group=tmesh)

##########################################################################################################
## Spatial index of the field component of the SPDE model
idx <- inla.spde.make.index('s', spde$n.spde, n.group=k)

##########################################################################################################
## Stack for fitting the model
stk <- inla.stack(data=list(y=agg.dat$Freq, exposure=e0),
                  A=list(A.st, 1),
                  effects=list(idx,
                               list(b0=rep(1, nrow(agg.dat)),rwtime=rwtime,seastime=seastime,artime=artime,
                                        vsolo=vsolo2,koppen=koppen,temp=vtemp,pluv=vpluv,declivsp=dummydecliv, 
                                        price=price, zoning=zone[,2])))

##########################################################################################################
## Fit model via INLA

## INLA formula
formula <- y ~ 0 + f(rwtime,model="rw1",constr=F)
                 +f(seastime,model="seasonal",season.length=4,constr=T)
                 +f(artime, constr = T, model = "ar", order = 2)
                 +as.factor(koppen)
                 +as.factor(zoning)
                 +as.factor(vsolo)
                 +price
                 +declivsp
                 +temp
                 +pluv
                 +f(s, model=spde, group=s.group, control.group=list(model='ar1'))

## Run INLA
res <- inla(formula, 
            family='poisson',
            data=inla.stack.data(stk),
            E=exposure,
            control.predictor=list(A=inla.stack.A(stk), compute=T),
            verbose=T,
            control.inla=list(correct=TRUE,correct.strategy='laplace'), 
            control.compute=list(config=T))


##########################################################################################################
## Trend component
library(Cairo)
library(xts)
trendm=as.xts(ts(res$summary.random$rwtime[,2],start = c(2003), frequency = 4))
trendl=as.xts(ts(res$summary.random$rwtime[,4],start = c(2003), frequency = 4))
trendu=as.xts(ts(res$summary.random$rwtime[,6],start = c(2003), frequency = 4))

trendmDF=data.frame(as.Date(index(trendm)),res$summary.random$rwtime[,2])
colnames(trendmDF)=c("Date","Trend_Mean")

trendlDF=data.frame(as.Date(index(trendm)),res$summary.random$rwtime[,4])
colnames(trendlDF)=c("Date","Trend_LI")

trenduDF=data.frame(as.Date(index(trendm)),res$summary.random$rwtime[,6])
colnames(trenduDF)=c("Date","Trend_LS")

trendciDF=data.frame(as.Date(index(trendm)),res$summary.random$rwtime[,4],res$summary.random$rwtime[,6])
colnames(trendciDF)=c("Date","a2.5q","a97.5q")

p1t <- ggplot() + 
  geom_line(data = trendmDF, aes(x = Date, y = Trend_Mean)) +
  xlab('') + scale_x_date(date_labels = "%Y",date_breaks="1 years")+   ylab('Trend')+
  geom_ribbon(data=trendciDF,aes(x=Date,ymin=a2.5q,ymax=a97.5q),alpha=0.3)

print(p1t)

##########################################################################################################
## Seasonal component
seasm=as.xts(ts(res$summary.random$seastime[,2],start = c(2003), frequency = 4))
seasl=as.xts(ts(res$summary.random$seastime[,4],start = c(2003), frequency = 4))
seasu=as.xts(ts(res$summary.random$seastime[,6],start = c(2003), frequency = 4))

seasmDF=data.frame(as.Date(index(seasm)),res$summary.random$seastime[,2])
colnames(seasmDF)=c("Date","seas_Mean")

seaslDF=data.frame(as.Date(index(seasm)),res$summary.random$seastime[,4])
colnames(seaslDF)=c("Date","seas_LI")

seasuDF=data.frame(as.Date(index(seasm)),res$summary.random$seastime[,6])
colnames(seasuDF)=c("Date","seas_LS")

seasciDF=data.frame(as.Date(index(seasm)),res$summary.random$seastime[,4],res$summary.random$seastime[,6])
colnames(seasciDF)=c("Date","a2.5q","a97.5q")

p1t <- ggplot() + 
  geom_line(data = seasmDF, aes(x = Date, y = seas_Mean)) +
  xlab('') + scale_x_date(date_labels = "%Y",date_breaks="1 years")+   ylab('seas')+
  geom_ribbon(data=seasciDF,aes(x=Date,ymin=a2.5q,ymax=a97.5q),alpha=0.3)

print(p1t)

##########################################################################################################
## Cyclic component
cyclem=as.xts(ts(res$summary.random$artime[,2],start = c(2003), frequency = 4))
cyclel=as.xts(ts(res$summary.random$artime[,4],start = c(2003), frequency = 4))
cycleu=as.xts(ts(res$summary.random$artime[,6],start = c(2003), frequency = 4))

cyclemDF=data.frame(as.Date(index(cyclem)),res$summary.random$artime[,2])
colnames(cyclemDF)=c("Date","cycle_Mean")

cyclelDF=data.frame(as.Date(index(cyclem)),res$summary.random$artime[,4])
colnames(cyclelDF)=c("Date","cycle_LI")

cycleuDF=data.frame(as.Date(index(cyclem)),res$summary.random$artime[,6])
colnames(cycleuDF)=c("Date","cycle_LS")

cycleciDF=data.frame(as.Date(index(cyclem)),res$summary.random$artime[,4],res$summary.random$artime[,6])
colnames(cycleciDF)=c("Date","a2.5q","a97.5q")

p1t <- ggplot() + 
  geom_line(data = cyclemDF, aes(x = Date, y = cycle_Mean)) +
  xlab('') + scale_x_date(date_labels = "%Y",date_breaks="1 years")+   ylab('cycle')+
  geom_ribbon(data=cycleciDF,aes(x=Date,ymin=a2.5q,ymax=a97.5q),alpha=0.3)

print(p1t)

## *LAURINI,  M.  A  spatio-temporal  approach  to  estimate  patterns  of  climate  change.
## Environmetrics, Wiley Online Library, v. 30, n. 1, p. e2542, 2019.
