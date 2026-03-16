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

##Data filter
#Install packages
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages<- c("dismo","terra","spatstat","cowplot","ggplot2","precrec","glmnet","maxnet","ranger")
check.packages(packages)
install.packages("mlr")


library(terra) #for spatial data
library(sf) #data frames with geometry
library(spatstat) #for point process modelling and converting between raster and pixel image objects

#Function to convert raster to images for spatstat. Takes one argument "im" that should be a raster object

raster.as.im = function(im) {
  #get the resolution (cell size of the raster)
  r = raster::res(im)[1]
  #get the origin (bottom left corner of the raster/image)
  orig = ext(im)[c(1,3)]
  #set the coordinates of the columns which is just a series of number from zero increasing by 100 metres (the resolution of the raster) for every cell along the rows and columns.
  xx = orig[1] + seq(from=0,to=(ncol(im) - 1)*100,by=r)
  #set the coordinates of the columns
  yy = orig[2] + seq(from=0,to=(nrow(im) - 1)*100,by=r)
  
  #now build a single matrix with the cell values and dimension we want - note that we reverse the rows in the matrix by setting nrow(im):1. This is just because spatstat reads images in reverse order to how raster object are usually organised (so the below code ensures our image is not upside-down when we come to analyse it)
  mat=matrix(raster::values(im), ncol = ncol(im), 
             nrow = nrow(im), byrow = TRUE)[nrow(im):1, ]
  return(spatstat.geom::im(mat, 
                           xcol = xx, yrow = yy))
}

#Use coordinates object to create our spatial points object
meles.sp=st_as_sf(meles.latlong,coords=c("x","y"),crs="epsg:4326")

#First set the extent to the study area
scot=st_read('scotSamp.shp')

#load in the land cover map and then clip to the polygon
LCM=rast("LCMUK.tif")

#crop to the extent of the study area plus a little more (because we will lose a small amount of data in the next step)
LCM=crop(LCM,st_buffer(scot, dist= 1000))

#aggregate LCM raster
LCM=aggregate(LCM$LCMUK_1,fact=4,fun="modal")
#project meles data
meles.sp=st_transform(meles.sp,crs(LCM))

#now crop our points to the study area
melesFin=meles.sp[scot,]

#finally, mask the LCM to this boundary
LCM=crop(LCM,scot,mask=TRUE)

#inspect
plot(LCM)
plot(melesFin$geometry,add=T)

#access levels of the raster by treating them as categorical data ('factors' in R)
LCM=as.factor(LCM$LCMUK_1)

#create an vector object called reclass
reclass = c(0,1,rep(0,20))

# combine with the LCM categories into a matrix of old and new values.
RCmatrix=cbind(levels(LCM)[[1]],reclass)
RCmatrix=RCmatrix[,2:3]

#apply function to make sure new columns are numeric (here the "2" specifies that we want to apply the as.numeric function to columns, where "1" would have specified rows)
RCmatrix=apply(RCmatrix,2,FUN=as.numeric)

#Use the classify() function to asssign new values to LCM with our reclassification matrix
broadleaf=classify(LCM, RCmatrix)

#neighbourhood weights matrix to sum all available resources for each cell
#get number of picels needed to cover the 900 metre radius (round() just rounds this to the nearest interger - because you can't have 0.5 pixels)

nPix=round(900/res(LCM)[1])

#next, you need to double this number (for the distance in two direction i.e. nPix is the radius and we want the diameter) and add one (because a weights matrix always needs to be odd so that you have a central cell and an equal number of grid cells either side).
nPix=(nPix*2)+1

#buiild weights matrix
weightsMatrix=matrix(1:nPix^2,nrow=nPix,ncol=nPix)

#get focal cell 
x=ceiling(ncol(weightsMatrix)/2)
y=ceiling(nrow(weightsMatrix)/2)
focalCell=weightsMatrix[x,y]
indFocal=which(weightsMatrix==focalCell,arr.ind = TRUE)

#compute distances
distances=list()

for(i in 1:nPix^2){
  ind.i=which(weightsMatrix==i,arr.ind=T)
  diffX=abs(ind.i[1,1]-indFocal[1,1])*res(LCM)[1]
  diffY=abs(ind.i[1,2]-indFocal[1,2])*res(LCM)[1]
  
  dist.i=sqrt(diffX^2+diffY^2)
  distances[[i]]=dist.i
  
}

#add distance values to the weights matrix
weightsMatrix[]=unlist(distances)

#set cells outside search radius to NA
weightsMatrix[weightsMatrix>900]=NA

#plot weights matrix
plot(rast(weightsMatrix))

#normalise the weights matrix by dividing all cell values by the number of cells. 
weightsMatrixNorm=weightsMatrix
weightsMatrixNorm[!is.na(weightsMatrixNorm)]=1/length(weightsMatrixNorm[!is.na(weightsMatrixNorm)])

#test to see for yourself
sum(weightsMatrixNorm,na.rm=T)

#plot weights matrix
plot(rast(weightsMatrixNorm))

#sum neighbourhood values from all surrounding cells
lcm_wood_900=focal(broadleaf,w=weightsMatrixNorm,fun="sum")
plot(lcm_wood_900)

#create an vector object called reclassUrban which is zero for all classes except tghe two urban classes in the LCM
reclassUrban = c(rep(0,19),1,1)

# combine with the LCM categories into a matrix of old and new values.
RCmatrixUrban= cbind(levels(LCM)[[1]],reclass)

RCmatrixUrban=RCmatrixUrban[,2:3]

#apply function to make sure new columns are numeric (here the "2" specifies that we want to apply the as.numeric function to columns, where "1" would have specified rows)
RCmatrixUrban=apply(RCmatrixUrban,2,FUN=as.numeric)
#Use the reclassify() function to asssign new values to LCM with our reclassification matrix
urban = classify(LCM, RCmatrixUrban)

#neighbourhood weights matrix to sum all available resources for each cell
#get number of picels needed to cover the 2300 metre radius for the urban class (round() just rounds this to the nearest interger - because you can't have 0.5 pixels)
nPixUrban=round(2300/res(LCM)[1])

#next, you need to double this number (one for rows and one for columns) and add one (because a weights matrix always needs to be odd so that you have a central cell and an equal number of grid cells either side).
nPixUrban=(nPixUrban*2)+1

#buiild weights matrix
weightsMatrixUrban=matrix(1:nPixUrban^2,nrow=nPixUrban,ncol=nPixUrban)

#get focal cell 
x=ceiling(ncol(weightsMatrixUrban)/2)
y=ceiling(nrow(weightsMatrixUrban)/2)
focalCell=weightsMatrixUrban[x,y]
indFocal=which(weightsMatrixUrban==focalCell,arr.ind = TRUE)

#compute distances
distancesUrban=list()

for(i in 1:nPixUrban^2){
  ind.i=which(weightsMatrixUrban==i,arr.ind=T)
  diffX=abs(ind.i[1,1]-indFocal[1,1])*res(LCM)[1]
  diffY=abs(ind.i[1,2]-indFocal[1,2])*res(LCM)[1]
  
  dist.i=sqrt(diffX^2+diffY^2)
  distancesUrban[[i]]=dist.i
  
}

#add distance values to the weights matrix
weightsMatrixUrban[]=unlist(distancesUrban)

#set cells outside search radius to NA
weightsMatrixUrban[weightsMatrixUrban>2300]=NA

#normalise the weights matrix by dividing all cell values by the number of cells. 
weightsMatrixUrban[!is.na(weightsMatrixUrban)]=1/length(weightsMatrixUrban[!is.na(weightsMatrixUrban)])

#sum urban class from all surrounding cells

lcm_urban_2300=focal(urban,w=weightsMatrixUrban,fun="sum")

#input the dem data
demScot=rast('demScotland.tif')
demScot=terra::resample(demScot,lcm_wood_900)

#inspect
plot(demScot)

#stack the covariate layers together
allEnv=c(lcm_wood_900,lcm_urban_2300,demScot)
names(allEnv)=c("broadleaf","urban","elev")

#creat background points
set.seed(11)

#sample background - one point for every cell (9775)
back = spatSample(allEnv,size=2000,as.points=TRUE,method="random",na.rm=TRUE) 
back=back[!is.na(back$broadleaf),]
back=st_as_sf(back,crs="EPSG:27700")

# get environmental covariates at presence locations
eP=terra::extract(allEnv,melesFin)

#bind together the presence data using cbind() which binds together objects by column
Pres.cov=st_as_sf(cbind(eP,melesFin))
Pres.cov$Pres=1

#Remove the first column which is just an ID field.
Pres.cov=Pres.cov[,-1]

#get coordinates for spatial cross-validation later
coordsPres=st_coordinates(Pres.cov)

#drop geometry column using st_drop_geometry()
Back.cov=st_as_sf(data.frame(back,Pres=0))


#get coordinates of background points for cross validation later
coordsBack=st_coordinates(back)

#combine
coords=data.frame(rbind(coordsPres,coordsBack))

#assign coumn names
colnames(coords)=c("x","y")

#combine pres and background
all.cov=rbind(Pres.cov,Back.cov)

#add coordinates
all.cov=cbind(all.cov,coords)

#remove any NAs
all.cov=na.omit(all.cov)
all.cov=st_drop_geometry(all.cov)

##Model evaluation
library(mlr)
#For the makeClassifTask function to work, our target variable needs to be categorical (a "factor" in R) so let's tidy that up first
task=all.cov
head(all.cov)
task$Pres=as.factor(task$Pres)
task = makeClassifTask(data = task[,c(1:4)], target = "Pres",
                       positive = "1", coordinates = task[,5:6])

## Binomial (logistic regression)
#use the make learner function to build the model approach. Fix factors prediction is set to TRUE here because our outcome is a factor (i.e. categorical: 0-1)
lrnBinomial = makeLearner("classif.binomial",
                          predict.type = "prob",
                          fix.factors.prediction = TRUE)

#set up resampling strategy (non-spatial cross-validation)
perf_levelCV = makeResampleDesc(method = "RepCV", predict = "test", folds = 5, reps = 5)

#set up resampling strategy for spatial cross-validation
perf_level_spCV = makeResampleDesc(method = "SpRepCV", folds = 5, reps = 5) #sampling strategy to run five fold re-sampling five times

#Binomial conventional cross validation (K fold)
cvBinomial = resample(learner = lrnBinomial, task =task,
                      resampling = perf_levelCV, 
                      measures = mlr::auc,
                      show.info = FALSE)
print(cvBinomial)

#create Spatial Resampling Plots
plots = createSpatialResamplingPlots(task,resample=cvBinomial,
                                     crs=crs(allEnv),datum=crs(allEnv),color.test = "red",point.size = 1)
library(cowplot)

#use the cowplot function to plot all folds out in a grid
cowplot::plot_grid(plotlist = plots[["Plots"]], ncol = 3, nrow = 2,
                   labels = plots[["Labels"]])

##Binomial spatial cross validation
sp_cvBinomial = resample(learner = lrnBinomial, task =task,
                         resampling = perf_level_spCV, 
                         measures = mlr::auc,
                         show.info = FALSE)
print(sp_cvBinomial)

#make partition plots
plotsSP = createSpatialResamplingPlots(task,resample=sp_cvBinomial,
                                       crs=crs(allEnv),datum=crs(allEnv),color.test = "red",point.size = 1)

#use the cowplot function to plot all folds out in a grid
cowplot::plot_grid(plotlist = plotsSP[["Plots"]], ncol = 3, nrow = 2,
                   labels = plotsSP[["Labels"]])

##Random Forest
lrnRF = makeLearner("classif.ranger",
                    predict.type = "prob",
                    fix.factors.prediction = TRUE)

#random sampling cross-validation
cvRF = resample(learner = lrnRF, task =task,
                resampling = perf_levelCV, 
                measures = mlr::auc,
                show.info = FALSE)
print(cvRF)

#spatial partitioning cross-validation
sp_cvRF = resample(learner = lrnRF, task =task,
                   resampling = perf_level_spCV, 
                   measures = mlr::auc,
                   show.info = FALSE)
print(sp_cvRF)

##tune the parameters of RF to provide more robust predictions to new, unknown locations
getParamSet(lrnRF)
paramsRF = makeParamSet(
  makeIntegerParam("mtry",lower = 1,upper = 3),
  makeIntegerParam("min.node.size",lower = 1,upper = 20),
  makeIntegerParam("num.trees",lower = 100,upper = 500)
)

# specifying random parameter value search
tune_level = makeResampleDesc(method = "SpCV", iters = 5)
ctrl = makeTuneControlRandom(maxit = 50)
tuned_RF = tuneParams(learner = lrnRF,
                      task = task,
                      resampling = tune_level,
                      measures = mlr::auc,
                      par.set = paramsRF,
                      control = ctrl,
                      show.info = FALSE)
print(tuned_RF)
