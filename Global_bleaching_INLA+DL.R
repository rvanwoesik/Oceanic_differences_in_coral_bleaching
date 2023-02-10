###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###-------------------- script for analyses of the global coral bleaching database --------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

## This code accompanies the paper titled "Oceanic differences in coral-bleaching responses to marine heatwaves"
## By Tom Shlesinger and Rob van Woesik
## Science of the Total Environment, 2023


## First, install and/or load the following packages:

library(tidyverse)
library(rworldmap)
library(colorspace)
library(sf)
library(INLA)
library(INLAutils)
library(ggregplot)
library(fields)
library(rgdal)
library(raster)
library(RColorBrewer)
library(magrittr)


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###--------------------------------------- Reading in the data ----------------------------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

setwd("") #set your working directory

home <- getwd()
 
data <- read.csv(file=file.path(home, "Bleaching_dataset.csv"), header=TRUE, sep=",")


#################################################################

# Choose the ocean to analyze 

#################################################################

#  Global

islsdat=data

# Define locations
Loc <- cbind(islsdat$Longitude_Degrees, islsdat$Latitude_Degrees)

#################################################################

#  Atlantic

islsdat=subset(data, Ocean_Name=="Atlantic" & Longitude_Degrees < 0) 

# Define locations
Loc <- cbind(islsdat$Longitude_Degrees, islsdat$Latitude_Degrees)


#################################################################

#  Pacific 

islsdat=subset(data, Ocean_Name=="Pacific") 

# Recentering to get the Pacific central:
# Convert data frame to sf object
my.sf.point <- st_as_sf(x = islsdat, 
                        coords = c("Longitude_Degrees", "Latitude_Degrees"),
                        crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

Loc1=st_shift_longitude(my.sf.point)
Loc=st_coordinates(Loc1)

#################################################################

#  Indian

islsdat=subset(data, Ocean_Name=="Indian") 

# Define locations
Loc <- cbind(islsdat$Longitude_Degrees, islsdat$Latitude_Degrees)



######################################################################
#Plot of the study sites 
wholeworld<-getMap(resolution="high")
land <- lighten("beige", amount = 0.4, space = "HCL")

plot(Loc, col = "red", xlab="",ylab="", cex.lab=1.5, cex.axis=1.2, type = "n")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "aliceblue")
plot(wholeworld, col = land, add = T)
par(new = T)
plot(Loc, col = "red", xlab="Longitude",ylab="Latitude", pch = 16, cex = 0.7, cex.lab=1.5, cex.axis=1.2)


######################################################################

# Standardize covariates in islsdat

standardize_function<-function(x){
  x.standardized=(x-mean(na.omit(x)))/sd(na.omit(x))
  return(x.standardized)
}

islsdat$SSTA_DHW<-standardize_function(islsdat$SSTA_DHW)
islsdat$SST<-standardize_function(islsdat$Temp)
islsdat$SSTA_Freq_stdev<-standardize_function(islsdat$SSTA_Frequency_Standard_Deviation)
islsdat$SST_stdev<-standardize_function(islsdat$Temperature_Kelvin_Standard_Deviation)
islsdat$TSA_DHW_stdev<-standardize_function(islsdat$TSA_DHW_Standard_Deviation)
islsdat$Turb <-standardize_function(islsdat$Turbidity)
islsdat$Dist <-standardize_function(islsdat$Distance_to_Shore)
islsdat$Depth<-standardize_function(islsdat$Depth_m)
islsdat$Coral_cover=standardize_function(islsdat$Coral)


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###---------------------------------------------- INLA ------------------------------------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

# Setting up mesh:

boundary <- inla.nonconvex.hull(points = Loc)
MeshC <- inla.mesh.2d(boundary = boundary, max.edge = c(2,6))
plot(MeshC, main="")
points(Loc, col = "red", pch = 16,  cex = 0.6)
MeshC$n

#After the mesh has been set up convert to model format.
#This uses an A matrix, which translates spatial locations on the mesh into vectors in the model.

A2 <- inla.spde.make.A(MeshC, loc = Loc)
dim(A2)  

###################################################
# Create SPDE
# The parameters "prior.range" and "prior.sigma" control the joint prior on range and standard deviation of the spatial field.
spde <- inla.spde2.pcmatern(mesh = MeshC, 
                                prior.range = c(50, 0.9), 
                                prior.sigma = c(1, 0.01))

#define spatial random field
w.index <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde,
  n.group = 1,
  n.repl  = 1)


##########################################################

# Create a data frame with an intercept and covariates.

N <- nrow(islsdat) 

X <- data.frame(Intercept = rep(1, N),
                Ocean = islsdat$Ocean_Name,
                Realm = islsdat$Realm_Name,
                SSTA_DHW  = islsdat$SSTA_DHW, 
                Site_ID=islsdat$Site_ID,
                SST=islsdat$SST,
                SSTA_Freq_stdev=islsdat$SSTA_Freq_stdev,
                SST_stdev=islsdat$SST_stdev,
                TSA_DHW_stdev=islsdat$TSA_DHW_stdev,
                Longitude=islsdat$Longitude_Degrees,
                Latitude=islsdat$Latitude_Degrees,
                Year=islsdat$Date_Year,
                Turb=islsdat$Turb,
                Dist=islsdat$Dist,
                Depth=islsdat$Depth,
                Coral_cover=islsdat$Coral_cover
)

##################################################################################

#                   Run spatial and temporal INLA Model

#################################################################################
# Full model with Beta distribution and fixed effects, random effects, and spatial random effects (SPDE)
# Spatial random effects taking care of temporal autocorrelation effects with a random-walk process

prec.prior <- list(prec = list(param = c(0.001, 0.001))) # for the hyper prior

# Smithson and Verkuilen 2006 transformation
islsdat$B=islsdat$avg.bleach/100
islsdat$B_1= (islsdat$B *(N-1)+0.5)/N

Stk10 <- inla.stack(
  tag  = "islsdat",
  data = list(y = islsdat$B_1),  
  A    = list(A2, 1),                      
  effects = list(                 
    w = w.index,            #Spatial field  
    X = as.data.frame(X)))  #Covariates


# For individual ocean analyses analysis use this formula:
F1 <- y ~ 1 + Coral_cover + Turb + Dist + Depth + SST + SSTA_DHW + SSTA_Freq_stdev+ SST_stdev + TSA_DHW_stdev + Longitude + Latitude +
  f(Year, model="rw1", hyper=prec.prior) + f(Site_ID, model="iid") + f(Realm, model="iid") +
  f(w, model=spde)


# For the global analysis use this formula:
F1 <- y ~ 1 + Coral_cover + Turb + Dist + Depth + SST + SSTA_DHW + SSTA_Freq_stdev+ SST_stdev + TSA_DHW_stdev + Longitude + Latitude +
  f(Year, model="rw1", hyper=prec.prior) + f(Site_ID, model="iid") + f(Ocean, model="iid") + f(Realm, model="iid") +
  f(w, model=spde)


## Run INLA:
I1 <- inla(F1,family = "beta",
            data = inla.stack.data(Stk10), verbose = T,
            control.compute = list(dic = TRUE, waic=TRUE, cpo=TRUE),
            control.predictor = list(A = inla.stack.A(Stk10), compute=TRUE))


summary(I1)
Efxplot(list(I1), Intercept = F) + theme_bw() + 
  theme(legend.position="top", panel.grid=element_blank())



#Validation 
autoplot(I1)
Fitted_values2= inla.stack.index(Stk10, tag = "islsdat")$data
Fit=I1$summary.fitted.values[Fitted_values2, "mean"]
plot(Fit, islsdat$B_1, xlab="Fitted values", ylab="Observed values", xlim=c(0,1))
abline(coef = c(0,1))
cor.test(Fit, islsdat$B_1)
hist(I1$cpo$pit, main="", breaks=100)


# Plot spatial effect + SE

Q = inla.spde.precision(spde, theta = c(0,0))
sample = inla.qsample(n = 2, Q)
proj <- inla.mesh.projector(MeshC, dims = c(300, 300))
sample_proj = inla.mesh.project(proj, field = I1$summary.random$w$mean)
wholeworld<-getMap(resolution="high")
land <- lighten("beige", amount = 0.4, space = "HCL")


## For the Indian Ocean:
image.plot(proj$x, proj$y, sample_proj, legend.args=list(text="Posteriori spatial latent effect",  
                                                         cex=1.2, side=4, line=2.5),legend.cex=1.2,
           xlim=c(20,128), ylim=c(-33,33),
           col=(tim.colors()),legend.shrink=1, xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land,add=TRUE); box()

sample_proj.s = inla.mesh.project(proj, field = I1$summary.random$w$sd)
image.plot(proj$x, proj$y, sample_proj.s, legend.args=list(text="Posteriori spatial latent effect SE", 
                                                           cex=1.2, side=4, line=2.5),legend.cex=1.2,
           xlim=c(20,128), ylim=c(-33,33),
           col=(tim.colors()),legend.shrink=1, xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land,add=TRUE); box()


## For the Atlantic:
image.plot(proj$x, proj$y, sample_proj, legend.args=list(text="Posteriori spatial latent effect",  
                                                         cex=1.2, side=4, line=2.5),legend.cex=1.2,
           xlim=c(-100,-30), ylim=c(-25,40),
           col=(tim.colors()),legend.shrink=1, xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land,add=TRUE); box()

sample_proj.s = inla.mesh.project(proj, field = I1$summary.random$w$sd)
image.plot(proj$x, proj$y, sample_proj.s, legend.args=list(text="Posteriori spatial latent effect SE", 
                                                           cex=1.2, side=4, line=2.5),legend.cex=1.2,
           xlim=c(-100,-30), ylim=c(-25,40),
           col=(tim.colors()),legend.shrink=1, xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land,add=TRUE); box()

## For the Pacific/Global:
image.plot(proj$x, proj$y, sample_proj, legend.args=list(text="Posteriori spatial latent effect",  
                                                         cex=1.2, side=4, line=2.5),legend.cex=1.2,
           col=(tim.colors()),legend.shrink=1, xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land,add=TRUE); box()

sample_proj.s = inla.mesh.project(proj, field = I1$summary.random$w$sd)
image.plot(proj$x, proj$y, sample_proj.s, legend.args=list(text="Posteriori spatial latent effect SE", 
                                                           cex=1.2, side=4, line=2.5),legend.cex=1.2,
           col=(tim.colors()),legend.shrink=1, xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land,add=TRUE); box()


##########################################################################

                   #Spatial predictions for coral bleaching

##########################################################################
#First generate a grid
ch <- chull(Loc)
coords <- Loc[c(ch, ch[1]), ]  # closed polygon

sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)), proj4string=CRS("+proj=longlat +datum=WGS84"))
SF <- SpatialPolygonsDataFrame(sp_poly, data=data.frame(ID=1))

writeOGR(SF, "hull", layer="SP", driver="ESRI Shapefile")
plot(Loc)
plot(SF, add=T)
path <-  paste(home, "/Hull", sep = "")
Carb <- readOGR(path, layer = "SP")
extent(Carb)
r <- raster(ncols=750, nrows=750)
CR <- rasterize(Carb, r)
y.res=CR@nrows
x.res=CR@ncols

#Create grid of ncol X nrow where we wish to project the model predictions 
Seq.X.grid= seq(from=CR@extent@xmin,
                to=CR@extent@xmax,
                length=x.res)
Seq.Y.grid= seq(from=CR@extent@ymin,
                to=CR@extent@ymax,
                length=y.res)
pred.grid=as.matrix(expand.grid(x=Seq.X.grid,
                                y=Seq.Y.grid))

####################################################
#Need 2 projection matrices and two stacks
MeshPred <- MeshC
spde.pred <- inla.spde2.matern(MeshPred, alpha = 2)

s.index.p<- inla.spde.make.index(name="sp.field.pred", n.spde=spde.pred$n.spde)

####################################################
A_est<- inla.spde.make.A(MeshPred, loc=Loc)
A_pred<- inla.spde.make.A(mesh=MeshPred)

Stk2 <- inla.stack(
  tag  = "Est",
  data = list(y = islsdat$B_1),  
  A    = list(A_est),                   
  effects = list(c(s.index.p, list(Intercept=1))))                 

Stk22 <- inla.stack(
  tag  = "Pred",
  data = list(y = NA),  
  A    = list(A_pred),                      
  effects = list(c(s.index.p, list(Intercept=1))))                 
  
StackJoin22=inla.stack(Stk2,Stk22)  

#############################################
#Inla

formula_Pred<- y ~ 1 + f(sp.field.pred, model=spde.pred)

Mod_Pred<- inla(formula_Pred, data = inla.stack.data(StackJoin22, spde = spde.pred),
                              family = "beta", verbose = T,
                              control.predictor = list(A = inla.stack.A(StackJoin22), compute = T))
summary(Mod_Pred)

###################################################################
# We need to extract the index of the data from the prediction part of the stack 
# (using the tag "Pred" we assigned to the stack) and use it to select the relevant 
# posterior mean and sd for the predicted response variable. 
# Then we use the inla.mesh.projector() function to calculate the projection from the 
# Mesh to the grid we created (pred.grid).

index.pred <- inla.stack.index(StackJoin22, "Pred")$data
post.mean.pred <- Mod_Pred$summary.linear.predictor[index.pred, "mean"]
post.sd.pred <- Mod_Pred$summary.linear.predictor[index.pred, "sd"]

#Plot on grid
proj.grid <- inla.mesh.projector(MeshPred,
                                 xlim = range(pred.grid[,1]),
                                 ylim = range(pred.grid[,2]), 
                                 dims = c(x.res, y.res))
                    
#Project the values we extracted from the model on the lattice we have created and transform the projected predictions to a raster object as we did before with the GRF and plot them in a similar fashion (we do this for both the mean and standard deviation).
post.mean.pred.grid <- inla.mesh.project(proj.grid, post.mean.pred)
post.sd.pred.grid <- inla.mesh.project(proj.grid, post.sd.pred)

predmean <- t(post.mean.pred.grid)
predmean2 <- predmean[rev(1:length(predmean[,1])),]

predmean_ras <- raster(predmean2,
                       xmn = range(proj.grid$x)[1], xmx = range(proj.grid$x)[2],
                       ymn = range(proj.grid$y)[1], ymx = range(proj.grid$y)[2])

predsd <- t(post.sd.pred.grid)
predsd2 <- predsd[rev(1:length(predsd[,1])),]
predsd_ras <- raster(predsd2,
                     xmn = range(proj.grid$x)[1], xmx = range(proj.grid$x)[2],
                     ymn = range(proj.grid$y)[1], ymx = range(proj.grid$y)[2])


# plot the model predictions for mean + sd
my.palette<-brewer.pal(n=9, name="YlOrRd")

# For the Indian Ocean:
plot(predmean_ras, asp = 0, col = my.palette, xlim=c(29,75), ylim=c(-32,6),
     xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land, add=TRUE,); box() ; axis(2, at = c(-30, -20, -10, 0, 10, 20, 30), cex.lab=1.5, cex.axis=1.2)

plot(predsd_ras, asp = 0, col = my.palette, xlim=c(29,75), ylim=c(-32,6),
     xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land,add=TRUE); box() ; axis(2, at = c(-30, -20, -10, 0, 10, 20, 30), cex.lab=1.5, cex.axis=1.2)


# For the Atlantic:
plot(predmean_ras, asp = 0, col = my.palette, xlim=c(-105,-27), ylim=c(-24,39),
     xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land, add=TRUE,); box() ; axis(2, at = c(-30, -20, -10, 0, 10, 20, 30), cex.lab=1.5, cex.axis=1.2)

plot(predsd_ras, asp = 0, col = my.palette, xlim=c(-105,-27), ylim=c(-24,39),
     xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land,add=TRUE); box() ; axis(2, at = c(-30, -20, -10, 0, 10, 20, 30), cex.lab=1.5, cex.axis=1.2)

# For the Pacific/Global:
plot(predmean_ras, asp = 0, col = my.palette, yaxt = "n", ylim=c(-40,40),
     xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land, add=TRUE,); box() ; axis(2, at = c(-30, -20, -10, 0, 10, 20, 30), cex.lab=1.5, cex.axis=1.2)

plot(predsd_ras, asp = 0, col = my.palette, yaxt = "n", ylim=c(-40,40),
     xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land,add=TRUE); box() ; axis(2, at = c(-30, -20, -10, 0, 10, 20, 30), cex.lab=1.5, cex.axis=1.2)

##################################################################


#    Spatial predictions for coral cover


##################################################################

# Smithson and Verkuilen 2006 transformation
islsdat$CC=islsdat$Coral/100
N=sum(!is.na(islsdat$CC))
islsdat$CC_1= (islsdat$CC *(N-1)+0.5)/N

#Need 2 projection matrices and two stacks
MeshPred <- MeshC
spde.pred <- inla.spde2.matern(MeshPred, alpha = 2)

s.index.p<- inla.spde.make.index(name="sp.field.pred", n.spde=spde.pred$n.spde)

####################################################
A_est<- inla.spde.make.A(MeshPred, loc=Loc)
A_pred<- inla.spde.make.A(mesh=MeshPred)

Stk2 <- inla.stack(
  tag  = "Est",
  data = list(y = islsdat$CC_1),  
  A    = list(A_est),                 
  effects = list(c(s.index.p, list(Intercept=1))))                 

Stk22 <- inla.stack(
  tag  = "Pred",
  data = list(y = NA),  
  A    = list(A_pred),                      
  effects = list(c(s.index.p, list(Intercept=1))))                 

StackJoin22=inla.stack(Stk2,Stk22)  

#############################################
#Inla

formula_Pred <- y ~ 1 + f(sp.field.pred, model=spde.pred)

Mod_Pred<- inla(formula_Pred, data = inla.stack.data(StackJoin22, spde = spde.pred),
                family = "beta", verbose = T,
                control.predictor = list(A = inla.stack.A(StackJoin22), compute = T))
summary(Mod_Pred)


###################################################################

#We need to extract the index of the data from the prediction part of the stack (using the tag "Pred" we assigned to the stack) and use it to select the relevant posterior mean and sd for the predicted response variable. Then we use the inla.mesh.projector() function to calculate the projection from the Mesh to the grid we created (pred.grid).
index.pred <- inla.stack.index(StackJoin22, "Pred")$data

post.mean.pred <- Mod_Pred$summary.linear.predictor[index.pred, "mean"]
post.sd.pred <- Mod_Pred$summary.linear.predictor[index.pred, "sd"]

#Plot on grid
proj.grid <- inla.mesh.projector(MeshPred,
                                 xlim = range(pred.grid[,1]),
                                 ylim = range(pred.grid[,2]), 
                                 dims = c(x.res, y.res))

#Project the values we extracted from the model on the lattice we have created and transform the projected predictions to a raster object as we did before with the GRF and plot them in a similar fashion (we do this for both the mean and standard deviation).
post.mean.pred.grid <- inla.mesh.project(proj.grid, post.mean.pred)
post.sd.pred.grid <- inla.mesh.project(proj.grid, post.sd.pred)

predmean <- t(post.mean.pred.grid)
predmean2 <- predmean[rev(1:length(predmean[,1])),]

predmean_ras <- raster(predmean2,
                       xmn = range(proj.grid$x)[1], xmx = range(proj.grid$x)[2],
                       ymn = range(proj.grid$y)[1], ymx = range(proj.grid$y)[2])

predsd <- t(post.sd.pred.grid)
predsd2 <- predsd[rev(1:length(predsd[,1])),]
predsd_ras <- raster(predsd2,
                     xmn = range(proj.grid$x)[1], xmx = range(proj.grid$x)[2],
                     ymn = range(proj.grid$y)[1], ymx = range(proj.grid$y)[2])


# plot the model predictions for mean + sd
my.palette<-brewer.pal(n=9, name="BuGn")


# For the Indian Ocean:

plot(predmean_ras, asp = 0, col = my.palette, xlim=c(29,75), ylim=c(-32,6),
     xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land, add=TRUE); box(); axis(2, at = c(-30, -15, 0, 15, 30), cex.lab=1.5, cex.axis=1.2)

plot(predsd_ras, asp = 0, col = my.palette, xlim=c(20,128), ylim=c(-33,33),
     xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land, add=TRUE); box(); axis(2, at = c(-30, -15, 0, 15, 30), cex.lab=1.5, cex.axis=1.2)


# For the Atlantic:

plot(predmean_ras, asp = 0, col = my.palette, xlim=c(-105,-27), ylim=c(-24,39),
     xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land, add=TRUE); box(); axis(2, at = c(-30, -15, 0, 15, 30), cex.lab=1.5, cex.axis=1.2)

plot(predsd_ras, asp = 0, col = my.palette, xlim=c(-105,-27), ylim=c(-24,39),
     xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2)
plot(wholeworld, col=land, add=TRUE); box(); axis(2, at = c(-30, -15, 0, 15, 30), cex.lab=1.5, cex.axis=1.2)

# For the Pacific/Global:

plot(predmean_ras, asp = 0, col = my.palette, yaxt = "n", 
     xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2, ylim=c(-40,40))
plot(wholeworld, col=land, add=TRUE); box(); axis(2, at = c(-30, -15, 0, 15, 30), cex.lab=1.5, cex.axis=1.2)

plot(predsd_ras, asp = 0, col = my.palette, yaxt = "n", 
     xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.axis=1.2,ylim=c(-40,40))
plot(wholeworld, col=land, add=TRUE); box(); axis(2, at = c(-30, -15, 0, 15, 30), cex.lab=1.5, cex.axis=1.2)





###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###----------------------------------- Machine Learning using H2O -------------------------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

#Set up H2O
if(!require(h2o)) install.packages("h2o")

h2o.init(nthreads = -1, #Number of threads -1 means use all cores on your machine
         max_mem_size = "8G")  #max mem size is the maximum memory to allocate to H2O

# Simplify variable names:
islsdat = data

islsdat$Turb <- islsdat$Turbidity
islsdat$Dist <- islsdat$Distance_to_Shore
islsdat$Depth <- islsdat$Depth_m
islsdat$SST <- islsdat$Temp
islsdat$SSTA_Freq_stdev <- islsdat$SSTA_Frequency_Standard_Deviation
islsdat$SST_stdev <- islsdat$Temperature_Kelvin_Standard_Deviation
islsdat$TSA_DHW_stdev <- islsdat$TSA_DHW_Standard_Deviation


##############################################################
#Extracting, preparing, and transferring the data from R to h2o

## For the Indian Ocean:
islsdat_ind=subset(islsdat, Ocean_Name=="Indian") 
data_ind<-as.h2o(islsdat_ind)

## For the Atlantic Ocean:
islsdat_atl=subset(islsdat, Ocean_Name=="Atlantic" & Longitude_Degrees < 0) 
data_atl<-as.h2o(islsdat_atl)

## For the Pacific Ocean:
islsdat_pac=subset(islsdat, Ocean_Name=="Pacific") 
data_pac<-as.h2o(islsdat_pac)

## For the global analysis:
data_glbl<-as.h2o(islsdat)


############################################################
# Define variables for the model

y <- "avg.bleach"
x <- c("TSA_DHW_stdev", "Turb", "Latitude_Degrees", "Longitude_Degrees", "Dist", "Coral", "SSTA_Freq_stdev", 
       "Depth", "SST", "SST_stdev", "SSTA_DHW")


############################################################
# Partition the data into training, validation and test sets
splits_atl <- h2o.splitFrame(data = data_atl, 
                             ratios = c(0.7, 0.15),  #partition data into 70%, 15%, 15% chunks
                             seed = 1)  #setting a seed will guarantee reproducibility
train_atl <- splits_atl[[1]]
valid_atl <- splits_atl[[2]]
test_atl <- splits_atl[[3]]

splits_ind <- h2o.splitFrame(data = data_ind, 
                             ratios = c(0.7, 0.15),  #partition data into 70%, 15%, 15% chunks
                             seed = 1)  #setting a seed will guarantee reproducibility
train_ind <- splits_ind[[1]]
valid_ind <- splits_ind[[2]]
test_ind <- splits_ind[[3]]

splits_pac <- h2o.splitFrame(data = data_pac, 
                             ratios = c(0.7, 0.15),  #partition data into 70%, 15%, 15% chunks
                             seed = 1)  #setting a seed will guarantee reproducibility
train_pac <- splits_pac[[1]]
valid_pac <- splits_pac[[2]]
test_pac <- splits_pac[[3]]

splits_glbl <- h2o.splitFrame(data = data_glbl, 
                              ratios = c(0.7, 0.15),  #partition data into 70%, 15%, 15% chunks
                              seed = 1)  #setting a seed will guarantee reproducibility
train_glbl <- splits_glbl[[1]]
valid_glbl <- splits_glbl[[2]]
test_glbl <- splits_glbl[[3]]

##############################################################################

Atlantic_model <- h2o.deeplearning(x = x,
                                   y = y,
                                   training_frame = train_atl,
                                   model_id = "Atlantic_model",
                                   validation_frame = valid_atl,  #in DL, early stopping is on by default
                                   epochs = 20,
                                   hidden = c(10,10),
                                   score_interval = 1,           #used for early stopping
                                   stopping_rounds = 3,          #used for early stopping
                                   stopping_metric = "deviance", #used for early stopping
                                   stopping_tolerance = 0.0005,  #used for early stopping
                                   seed = 1,
                                   adaptive_rate = TRUE,
                                   input_dropout_ratio = 0.1,
                                   overwrite_with_best_model = TRUE)

Atlantic_model

dl_perf3 <- h2o.performance(model = Atlantic_model,
                            newdata = test_atl)
dl_perf3

par(mfrow=c(3,4))
x_atl <- h2o.partialPlot(object = Atlantic_model, data = data_atl, 
                         cols = c("TSA_DHW_stdev", "Turb", "Latitude_Degrees", "Longitude_Degrees", "Dist", "Coral", "SSTA_Freq_stdev", "Depth", "SST", "SST_stdev", "SSTA_DHW"))


Indian_model <- h2o.deeplearning(x = x,
                                 y = y,
                                 training_frame = train_ind,
                                 model_id = "Indian_model",
                                 validation_frame = valid_ind,  #in DL, early stopping is on by default
                                 epochs = 20,
                                 hidden = c(10,10),
                                 score_interval = 1,           #used for early stopping
                                 stopping_rounds = 3,          #used for early stopping
                                 stopping_metric = "deviance", #used for early stopping
                                 stopping_tolerance = 0.0005,  #used for early stopping
                                 seed = 1,
                                 adaptive_rate = TRUE,
                                 input_dropout_ratio = 0.1,
                                 overwrite_with_best_model = TRUE)

Indian_model

dl_perf3 <- h2o.performance(model = Indian_model,
                            newdata = test_ind)
dl_perf3

x_ind <- h2o.partialPlot(object = Indian_model, data = data_ind,
                         cols = c("TSA_DHW_stdev", "Turb", "Latitude_Degrees", "Longitude_Degrees", "Dist", "Coral", "SSTA_Freq_stdev", "Depth", "SST", "SST_stdev", "SSTA_DHW"))


Pacific_model <- h2o.deeplearning(x = x,
                                  y = y,
                                  training_frame = train_pac,
                                  model_id = "Pacific_model",
                                  validation_frame = valid_pac,  #in DL, early stopping is on by default
                                  epochs = 20,
                                  hidden = c(10,10),
                                  score_interval = 1,           #used for early stopping
                                  stopping_rounds = 3,          #used for early stopping
                                  stopping_metric = "deviance", #used for early stopping
                                  stopping_tolerance = 0.0005,  #used for early stopping
                                  seed = 1,
                                  adaptive_rate = TRUE,
                                  input_dropout_ratio = 0.1,
                                  overwrite_with_best_model = TRUE)

Pacific_model

dl_perf3 <- h2o.performance(model = Pacific_model,
                            newdata = test_pac)
dl_perf3

x_pac <- h2o.partialPlot(object = Pacific_model, data = data_pac,
                         cols = c("TSA_DHW_stdev", "Turb", "Latitude_Degrees", "Longitude_Degrees", "Dist", "Coral", "SSTA_Freq_stdev", "Depth", "SST", "SST_stdev", "SSTA_DHW"))


Global_model <- h2o.deeplearning(x = x,
                                 y = y,
                                 training_frame = train_glbl,
                                 model_id = "Global_model",
                                 validation_frame = valid_glbl,  #in DL, early stopping is on by default
                                 epochs = 20,
                                 hidden = c(10,10),
                                 score_interval = 1,           #used for early stopping
                                 stopping_rounds = 3,          #used for early stopping
                                 stopping_metric = "deviance", #used for early stopping
                                 stopping_tolerance = 0.0005,  #used for early stopping
                                 seed = 1,
                                 adaptive_rate = TRUE,
                                 input_dropout_ratio = 0.1,
                                 overwrite_with_best_model = TRUE)

Global_model

dl_perf3 <- h2o.performance(model = Global_model,
                            newdata = test_glbl)
dl_perf3

x_glbl <- h2o.partialPlot(object = Global_model, data = data_glbl,
                          cols = c("TSA_DHW_stdev", "Turb", "Latitude_Degrees", "Longitude_Degrees", "Dist", "Coral", "SSTA_Freq_stdev", "Depth", "SST", "SST_stdev", "SSTA_DHW"))
