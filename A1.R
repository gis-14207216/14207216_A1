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

#create an vector object called reclass
reclass <- c(0,1,rep(0,19))

# combine with the LCM categories into a matrix of old and new values.
RCmatrix<- cbind(levels(LCM)[[1]],reclass)

RCmatrix<-RCmatrix[,2:3]

#apply function to make sure new columns are numeric (here the "2" specifies that we want to apply the as.numeric function to columns, where "1" would have specified rows)
RCmatrix=apply(RCmatrix,2,FUN=as.numeric)
#Use the reclassify() function to asssign new values to LCM with our reclassification matrix
RCmatrix

broadleaf <- classify(LCM, RCmatrix)

#plot
plot(broadleaf)
plot(melesFin,add=TRUE)

#function for automating whole dataset.

landBuffer <- function(speciesData, r){         
  
  #buffer each point
  melesBuffer <- st_buffer(speciesData, dist=r)                     
  
  #crop the woodland layer to the buffer extent
  bufferlandcover <- crop(broadleaf, melesBuffer)              
  
  # now extract the raster values (which should all be 1 for woodland and 0 for everything else) within each buffer and sum to get number of woodland cells inside the buffers.
  masklandcover <- extract(bufferlandcover, melesBuffer,fun="sum")      
  #get woodland area (625 is the area in metres of each cell of our 25m raster)
  landcoverArea <- masklandcover$LCMUK_1*625  
  
  # convert to precentage cover (we use the st_area() function from the sf package to get the area of our buffer) but convert to a numeric object (because sf applies units i.e. metres which then cant be entered into numeric calculations)
  percentcover <- landcoverArea/as.numeric(st_area(melesBuffer))*100 
  
  # return the result
  return(percentcover)                                       
}
resList=list()

#loop
radii<-seq(100,2000,by=100)
for(i in radii){
  res.i=landBuffer(speciesData=melesSF,r=i)
  res.i
  resList[[i/100]]=res.i
  print(i)
  
}

#collect all results together
resFin=do.call("cbind",resList)

#convert to data frame
glmData=data.frame(resFin)

#assign more intuitive column names
colnames(glmData)=paste("radius",radii,sep="")

head(glmData)

#add in the presences data
glmData$Pres<-melesData$Pres

head(glmData)

#init empty data frame
glmRes=data.frame(radius=NA,loglikelihood=NA)

#for loop to iterate over radius values and run a general linear model with glm()
for(i in radii){
  
  #build the model formula  
  n.i=paste0("Pres~","radius",i,sep ="")
  
  #run
  glm.i=glm(formula(n.i),family = "binomial",data = glmData)
  
  #get logliklihood
  ll.i=as.numeric(logLik(glm.i))
  
  #collect results
  glmRes=rbind(glmRes,c(i,ll.i))
  
}

#inspect
head(glmRes)

plot(glmRes$radius, glmRes$loglikelihood,
     type="b",
     pch=19,
     lwd=2,
     col="red",
     xlab="Buffer radius (m)",
     ylab="Log-likelihood",
     cex.lab=1.3,
     mgp=c(3.0,1,0))

box()

#determine the optimum buffer size

#remove the NAs in the first row
glmRes=glmRes[!is.na(glmRes),]

#use the which.max function to subset the dataframe to the just the row containing the max log likelihood value

opt=glmRes[which.max(glmRes$loglikelihood),]

#print
print(opt)