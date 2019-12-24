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

setwd("C:/Users/ferna/Documents/Doc USP/Espacial/sugarcane")
#setwd("f:/dados/modisfires")

fires <- readOGR("DL_FIRE_M6_19834/fire_archive_M6_19834.shp", layer="fire_archive_M6_19834") 
firesc=fires[fires$CONFIDENCE==100 & fires$SATELLITE=="Aqua" & substr(fires$ACQ_DATE,1,4)<2018,]
firesc=fires[fires$CONFIDENCE==100 & fires$SATELLITE=="Aqua",]
firesc=fires[fires$CONFIDENCE>50 & substr(fires$ACQ_DATE,1,4)<2017 & substr(fires$ACQ_DATE,1,4)>2002,]

data1=as.matrix(firesc$ACQ_DATE)
idata1=order(data1)
firesco=firesc[idata1,]

tsol2003 <- raster("2003_01_01.LC_Type1.tif") 
tsol2004 <- raster("2004_01_01.LC_Type1.tif") 
tsol2005 <- raster("2005_01_01.LC_Type1.tif") 
tsol2006 <- raster("2006_01_01.LC_Type1.tif") 
tsol2007 <- raster("2007_01_01.LC_Type1.tif") 
tsol2008 <- raster("2008_01_01.LC_Type1.tif") 
tsol2009 <- raster("2009_01_01.LC_Type1.tif") 
tsol2010 <- raster("2010_01_01.LC_Type1.tif") 
tsol2011 <- raster("2011_01_01.LC_Type1.tif") 
tsol2012 <- raster("2012_01_01.LC_Type1.tif") 
tsol2013 <- raster("2013_01_01.LC_Type1.tif") 
tsol2014 <- raster("2014_01_01.LC_Type1.tif") 
tsol2015 <- raster("2015_01_01.LC_Type1.tif") 
tsol2016 <- raster("2016_01_01.LC_Type1.tif") 

shapeest=readOGR("estadosl_2007.shp", layer="estadosl_2007") 
shapesp=shapeest[shapeest$SIGLAUF3=="SP",]

border=shapesp@polygons
bordersp=border[[1]]@Polygons[[1]]@coords

require(splancs)
xy.in <- inout(coordinates(firesco), bordersp)
firesco=firesco[xy.in,]

n<- length(firesco)

meses=as.matrix(substr(firesco$ACQ_DATE,1,7))
umeses=unique(meses)

quarter=meses
for (j in 1:length(meses)){
  if (substr(meses[j],6,7)=="01" | substr(meses[j],6,7)=="02" |substr(meses[j],6,7)=="03") quarter[j]=1 
  if (substr(meses[j],6,7)=="04" | substr(meses[j],6,7)=="05" |substr(meses[j],6,7)=="06") quarter[j]=2 
  if (substr(meses[j],6,7)=="07" | substr(meses[j],6,7)=="08" |substr(meses[j],6,7)=="09") quarter[j]=3 
  if (substr(meses[j],6,7)=="10" |substr(meses[j],6,7)=="11"|substr(meses[j],6,7)=="12") quarter[j]=4 
}

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

locs=coordinates(firesco)
bounds <- gUnaryUnion(shapesp)

smesh <- inla.mesh.2d(boundary=inla.sp2segment(bounds),
                      max.edge=1, cutoff=0.3)
library(Cairo)
Cairo(file="smesh.png",
      type = "png")
plot(smesh, frame.plot=F)
dev.off()

ndays= max(time)

k <- ndays; tmesh <- inla.mesh.1d(seq(0, ndays, length=k))

library(deldir)
dd <- deldir(smesh$loc[,1], smesh$loc[,2])
tiles <- tile.list(dd)

polys <- SpatialPolygons(lapply(1:length(tiles), function(i)
{ p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
n <- nrow(p)
Polygons(list(Polygon(p[c(1:n, 1),])), i)
}))

area <- factor(over(SpatialPoints(cbind(locs[,1], locs[,2])),
                    polys), levels=1:length(polys))

png("area.png")
plot(area)
dev.off()

#summary(summary(area))

otime=time
t.breaks <- sort(c(tmesh$loc[c(1,k)],tmesh$loc[2:k-1]/2 + tmesh$loc[2:k]/2))
table(time <- factor(findInterval(time, t.breaks),levels=1:(length(t.breaks)-1)))
time <- factor(otime,levels=1:(length(t.breaks)-1))

#summary(summary(time))

png("time.png")
plot(time)
dev.off()

agg.dat <- as.data.frame(table(area, time))

for(j in 1:2) ### set time and area as integer
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
#The space-time volumn (area unit per time unit) at each polygon and time knot
summary(e0 <- w.areas[agg.dat$area] * (w.t[agg.dat$time]))

A.st <- inla.spde.make.A(smesh, smesh$loc[agg.dat$area,],group=agg.dat$time, mesh.group=tmesh)
spde <- inla.spde2.matern(smesh)
idx <- inla.spde.make.index('s', spde$n.spde, n.group=k)

rwtime=agg.dat$time
seastime=agg.dat$time
artime=agg.dat$time

loci1=smesh$loc[agg.dat$area,1:2]
vtime=agg.dat$time

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

vsolo=c(vsolo2003,vsolo2004,vsolo2005,vsolo2006,vsolo2007,vsolo2008,vsolo2009,vsolo2010,vsolo2011,vsolo2012,vsolo2013,vsolo2014,vsolo2015,vsolo2016)
vsolo2=vsolo==12*1

flk=raster("Koppen_Brazil_2013/Koppen Brazil 2013/Raster/koppen_paper/w001000.adf")
koppen=extract(flk,loci1)

e <- extent(shapesp)
flc<- crop(flk, e)

Cairo(file="koppenimage.png", type="png")
plot(flc,legend=F,col = rev(terrain.colors(12)), frame.plot=F)
plot(bounds, add=T, frame.plot=F)
title('Koppen Climate Classification')
op <- par(cex = 0.8)
legend("bottomright", legend = c("Cwa","Am","Af","Cfa","Cwb","Csb","Csa","Cfb","BSh","As","Cwc","Aw"), fill = rev(terrain.colors(12)), adj=0.1)
dev.off()

vtemp=c()
for (jj in 1:56){
  tempbr= raster("rastersudestetempmax.tif",band=168+jj) 
  indextemp=vtime==jj
  ct1=loci1[indextemp,]
  vtempi=extract(tempbr,ct1)
  vtemp=c(vtemp,vtempi)
}

Cairo(file="temp2016_4.png", type="png")
plot(tempbr)
plot(bounds, add=T)
title('Max. Temperatures, Southeast Region, 2016Q4')
dev.off()

vpluv=c()
for (jj in 1:56){
  tempbr= raster("rastersudestepluv.tif",band=168+jj) 
  indextemp=vtime==jj
  ct1=loci1[indextemp,]
  vpluvi=extract(tempbr,ct1)
  vpluv=c(vpluv,vpluvi)
}

Cairo(file="pluv2016_4.png", type="png")
plot(tempbr)
plot(bounds, add=T)
title('Precipitation, Southeast Region, 2016Q4')
dev.off()

Cairo(file='modislandclass.png', type="png")
par(mfrow=c(2,2))
plot(tsol2004, main='Modis Land Cover Classification 2004')
plot(bounds, add=T)
plot(tsol2008, main='Modis Land Cover Classification 2008')
plot(bounds, add=T)
plot(tsol2012, main='Modis Land Cover Classification 2012')
plot(bounds, add=T)
plot(tsol2016, main='Modis Land Cover Classification 2016')
plot(bounds, add=T)
dev.off()

decliv <- raster("declividade_br.asc")
declivsp<-extract(decliv,loci1)
dummydecliv = declivsp
dummydecliv[]<- ifelse(declivsp[]>12,1,0)

e <- extent(bounds)
declsp<- crop(decliv, e)
dummydeclsp = declsp
dummydeclsp[]<-ifelse(declsp[]>12,1,0)

Cairo(file='declive.png', type="png")
plot(declsp)
plot(bounds, add=T)
title('Soil Slope')
dev.off()

zoneamento<- readOGR("ZoneamentoCana_4C.shp")
zone<-extract(zoneamento,loci1)

plotfire2016<-(firesco[firesco$CONFIDENCE>50 & substr(firesco$ACQ_DATE,1,4)==2016,])
plotfire2003<-(firesco[firesco$CONFIDENCE>50 & substr(firesco$ACQ_DATE,1,4)==2003,])
png('plot20162003.png')
par(mfrow=c(2,1))
plot(declsp, main="Soil Slope and Fires, 2016")
plot(bounds,add=T)
points(plotfire2016, pch=".")
plot(declsp, main="Soil Slope and Fires, 2003")
plot(bounds,add=T)
points(plotfire2003, pch=".")
dev.off()

png("declivsp.png")
plot(declsp,legend=T)
plot(bounds, add=T)
title('Soil Slope')
dev.off()

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

Cairo(file = "price.png", type = "png")
plot(ps)
dev.off()

A.st <- inla.spde.make.A(smesh, smesh$loc[agg.dat$area,],
                         group=agg.dat$time, mesh.group=tmesh)
spde <- inla.spde2.matern(smesh)
idx <- inla.spde.make.index('s', spde$n.spde, n.group=k)

stk <- inla.stack(data=list(y=agg.dat$Freq, exposure=e0),
                  A=list(A.st, 1),
                  effects=list(idx,list(b0=rep(1, nrow(agg.dat)),rwtime=rwtime,seastime=seastime,artime=artime,vsolo=vsolo2,koppen=koppen,temp=vtemp,pluv=vpluv,declivsp=dummydecliv, price=price, zoneamento=zone[,2])))


#1
formula <- y ~ 0 + f(rwtime,model="rw1",constr=F)+as.factor(vsolo)+f(seastime,model="seasonal",season.length=4,constr=T)+f(artime, constr = T, model = "ar", 
                                                                                                                           order = 2)+f(s, model=spde, group=s.group, control.group=list(model='ar1'))
#2
formula <- y ~ 0 + f(rwtime,model="rw2",constr=F)+as.factor(koppen)+as.factor(vsolo)+f(seastime,model="seasonal",season.length=4,constr=T)+f(artime, constr = T, model = "ar", 
                                                                                                                                             order = 2)+f(s, model=spde, group=s.group, control.group=list(model='ar1'))
#3
formula <- y ~ 0 + f(rwtime,model="rw2",constr=F)+declivsp+as.factor(koppen)+as.factor(vsolo)+f(seastime,model="seasonal",season.length=4,constr=T)+f(artime, constr = T, model = "ar", 
                                                                                                                                                      order = 2)+f(s, model=spde, group=s.group, control.group=list(model='ar1'))
#4
formula <- y ~ 0 + f(rwtime,model="rw1",constr=F)+as.factor(koppen)+as.factor(vsolo)+temp+pluv+f(seastime,model="seasonal",season.length=4,constr=T)+f(artime, constr = T, model = "ar", 
                                                                                                                                                       order = 2)+f(s, model=spde, group=s.group, control.group=list(model='ar1'))
#5
formula <- y ~ 0 + f(rwtime,model="rw1",constr=F)+as.factor(koppen)+as.factor(zoneamento)+as.factor(vsolo)+price+declivsp+temp+pluv+f(seastime,model="seasonal",season.length=4,constr=T)+f(artime, constr = T, model = "ar", 
                                                                                                                                                                      order = 2)+f(s, model=spde, group=s.group, control.group=list(model='ar1'))
#5
formula <- y ~ 0 + f(rwtime,model="rw1",constr=F)+as.factor(koppen)+as.factor(zoneamento)+as.factor(vsolo)+price+declivsp+temp+pluv+f(seastime,model="seasonal",season.length=4,constr=T)+f(artime, constr = T, model = "ar", 
                                                                                                                                                                                            order = 2)#+f(s, model=spde, group=s.group, control.group=list(model='ar1'))
#6
formula <- y ~ 0 + f(rwtime,model="ar1",constr=F)+f(seastime,model="seasonal",season.length=4,constr=T)+f(artime, constr = T, model = "ar", 
                                                                                                          order = 2)+f(s, model=spde, group=s.group, control.group=list(model='ar1'))
#7
formula <- y ~ 0 + f(rwtime,model="rw1",constr=F,hyper = list(theta1 = list(initial = -0.4509856, fixed = F)))+declivsp+as.factor(koppen)+as.factor(vsolo)+temp+pluv+f(seastime,model="seasonal",season.length=4,constr=T)+f(artime, constr = T, model = "ar", 
                                                                                                                                                                                                                             order = 2)+f(s, model=spde, group=s.group, control.group=list(model='ar1'))
#8
formula <- y ~ 0 + f(rwtime,model="rw1",constr=F,hyper = list(theta1 = list(initial = -0.4509856, fixed = F)))+as.factor(koppen)+as.factor(vsolo)+temp+pluv+f(seastime,model="seasonal",season.length=4,constr=T)+f(artime, constr = T, model = "ar", 
                                                                                                                                                                                                                    order = 2)+f(s, model=spde,  group=s.group, control.group=list(model='ar1'))
#9
formula <- y ~ 0 + f(rwtime,model="rw1",constr=F,hyper = list(theta1 = list(initial = -0.4509856, fixed = F)))+declivsp+as.factor(koppen)+as.factor(vsolo)+temp+pluv+f(seastime,model="seasonal",season.length=4,constr=T)+f(artime, constr = T, model = "ar", 
                                                                                                                                                                                                                             order = 2)

res <- inla(formula, family='poisson',
            data=inla.stack.data(stk), E=exposure,
            control.predictor=list(A=inla.stack.A(stk), compute=T),verbose=T,control.inla=list(correct=TRUE,correct.strategy='laplace'), control.compute=list(config=T))

summary(res)

library(excursions)

cb <- simconf.inla(res, stk, tag = "rwtime", link=TRUE, alpha=0.05)

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
Cairo(file = "trend.png", type="png")
print(p1t)
dev.off()

p1t <- ggplot() + geom_line(data = ps, aes(x = Date, y = Sugar_Price))
print(p1t)

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

Cairo(file = "seas.png", type="png")
print(p1t)
dev.off()

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

Cairo(file = "cycle.png", type="png")
print(p1t)
dev.off()

prj <- inla.mesh.projector(smesh, xlim=bbox(bounds)[1,],
                           ylim=bbox(bounds)[2,], dims=c(700, 700))
g.no.in <- is.na(over(SpatialPoints(prj$lattice$loc), bounds))

t.mean <- lapply(1:k, function(j) {
  z <- inla.mesh.project(prj, res$summary.ran$s$mean[idx$s.group==j])
  z[g.no.in] <- NA
  return(z)
})

fit.mean <- lapply(1:k, function(j) {
  z <- inla.mesh.project(prj, res$summary.ran$s$mean[idx$s.group==j])+res$summary.random$rwtime[j,2]+res$summary.random$seastime[j,2]+res$summary.random$artime[j,2]
  z[g.no.in] <- NA
  return(z)
})

zlims <- range(unlist(t.mean), na.rm=TRUE)
library(fields)
png('randomfieldcana1.png')
par(mfrow=c(8,4), mar=c(0,0,0,0))
for (j in 1:28) {
  image(prj$x, prj$y, t.mean[[j]],
        axes=FALSE, zlim=zlims, col=tim.colors(30))
  lines(bounds)
  #points(locs[time==j,1], locs[time==j,2], pch=".")
}
image.plot(prj$x, prj$y, t.mean[[j]]+1e9, axes=FALSE, zlim=zlims, xlab='',
           legend.mar=10, legend.width=5, col=tim.colors(30), horizontal=T)
dev.off()
##
zlims <- range(unlist(fit.mean), na.rm=TRUE)
library(fields)
png("fitcana2.png")
par(mfrow=c(8,4), mar=c(0,0.1,0,0.1))
for (j in 29:56) {
  image(prj$x, prj$y, fit.mean[[j]],
        axes=FALSE, zlim=zlims, col=tim.colors(30))
  lines(bounds)
  points(locs[time==j,1], locs[time==j,2], pch=".")
}  
image.plot(prj$x, prj$y, fit.mean[[j]]+1e9, axes=FALSE, zlim=zlims, xlab='',
           legend.mar=10, legend.width=5, col=tim.colors(30), horizontal=T)
dev.off()
library(Cairo)
Cairo(file = "zone.png", type="png")
plot(zoneamento,col = rev(terrain.colors(12)))
legend("bottomleft", legend = c("Dam","Appropriate","Appropriate - limitations","Appropriate - restrictions", "Inappropriate"), fill = rev(terrain.colors(12)))
dev.off()


roots= inla.ar.pacf2phi(c(res$summary.hyper[4,1],res$summary.hyper[5,1]))


(roots[1]^2+4*roots[2])

(roots[1]^2+4*roots[2])<0

f0=(1/(2*pi))*acos(roots[1]/(2*sqrt(-roots[2])))

1/f0




















lon.prd <- c(-53,-45)
lat.prd <- c(-26,-20)
stepsize <- 1
nxy <- round(c(diff(range(lon.prd)),diff(range(lat.prd)))/stepsize)

lattice.prd <- inla.mesh.projector(smesh,xlim=range(lon.prd),ylim=range(lat.prd), dims=nxy+1)
lattice.prd <- submesh.grid(matrix(1,nxy[1]+1,nxy[2]+1),list(loc=lattice.prd$lattice$loc,dims=nxy+1))

A.prd <- inla.spde.make.A(smesh, loc= loci1)

A.st <- inla.spde.make.A(smesh, smesh$loc[agg.dat$area,],group=agg.dat$time, mesh.group=tmesh)

spde <- inla.spde2.matern(smesh)

idx <- inla.spde.make.index('s', spde$n.spde, n.group=k)
mesh.index <- inla.spde.make.index(name = "field", n.spde = spde$n.spde)

stk.dat <- inla.stack(data=list(y=agg.dat$Freq, exposure=e0),
                  A=list(A.st, 1),
                  effects=list(idx,list(b0=rep(1, nrow(agg.dat)),rwtime=rwtime,seastime=seastime,artime=artime,vsolo=vsolo2,koppen=koppen,temp=vtemp,pluv=vpluv,declivsp=dummydecliv, price=price, zoneamento=zone[,2])))

stk.prd <- inla.stack(data=list(y = NA), A = list(A.prd, 1), tag = "prd", effects=list(c(mesh.index, list(rwtime=rep(1,smesh$n))),list(lat = submesh$loc[,2], lon = submesh$loc[,1])))

stk <- inla.stack(stk.dat, stk.prd)

formula <- y ~ 0 + f(rwtime,model="rw1",constr=F)+as.factor(koppen)+as.factor(zoneamento)+as.factor(vsolo)+price+declivsp+temp+pluv+f(seastime,model="seasonal",season.length=4,constr=T)+f(artime, constr = T, model = "ar", 
                                                                                                                                                                                            order = 2)+f(s, model=spde, group=s.group, control.group=list(model='ar1'))
res <- inla(formula, family='poisson',
            data=inla.stack.data(stk), E=exposure,
            control.predictor=list(A=inla.stack.A(stk), compute=T),verbose=T,control.inla=list(correct=TRUE,correct.strategy='laplace'), control.compute=list(config=T))




