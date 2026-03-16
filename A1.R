setwd("/Users/leo/Desktop/GEOG79122")

##Data input
#Install packages
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages<- c("terra","sf","mapview")
install.packages("terra")
check.packages(packages)

library(terra)
library(sf)

#read in the species occurence data
meles = read.csv("Melesmeles.csv")

#subset the data to only include points with complete coordinates.
meles=meles[!is.na(meles$Latitude),]

#remove all points with uncertainty > 1000m
meles<-meles[meles$Coordinate.uncertainty_m<=1000,]

#make spatial points layer
meles.latlong=data.frame(x=meles$Longitude,y=meles$Latitude)

#Use coordinates object to create our spatial points object
meles.sp<-vect(meles.latlong,geom=c("x","y"))

#check that the points now have our desired crs. 
crs(meles.sp)<-"epsg:4326"
plot(meles.sp)
#set the exetent to something workable
studyExtent=c(-4.2,-2.7,56.5,57.5) #list coordinates in the order: min x, max x, min y, max y

# crop points to study area
C=crop(meles.sp,studyExtent)

#read in the raster data
LCM=rast("LCMUK.tif")

#project points to land cover data
melesFin<-project(C,crs(LCM))

#get coordinates and set extent with 5km buffer
melesCoords<-crds(melesFin)

x.min <- min(melesCoords[,1]) - 5000
x.max <- max(melesCoords[,1]) + 5000
y.min <- min(melesCoords[,2]) - 5000
y.max <- max(melesCoords[,2]) + 5000
extent.new <- ext(x.min, x.max, y.min, y.max)

#crop the LCM to this extent
LCM <- crop(LCM$LCMUK_1, extent.new)

#check results
plot(LCM)
plot(melesFin,add=TRUE)

set.seed(11)

back.xy <- spatSample(LCM, size=1000,as.points=TRUE) 

#create a spatialPoints layer from the back.xy matrix

plot(LCM)
plot(melesFin,add=T)
plot(back.xy,add=TRUE, col='red')

#extract LCM data to presnec and background points
eA=extract(LCM,back.xy)
eP=extract(LCM,melesFin)
Abs=data.frame(crds(back.xy),Pres=0)
Pres=data.frame(crds(melesFin),Pres=1)

#check data
head(Pres)
head(Abs)
# bind the two data frames by row (both dataframes have the same column headings)
melesData=rbind(Pres,Abs)

#inspect
head(melesData)
melesSF=st_as_sf(melesData,coords=c("x","y"),crs="EPSG:27700")

#access levels of the raster by treating them as categorical data
LCM<-as.factor(LCM)
levels(LCM)