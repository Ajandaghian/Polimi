#_________________ Applied Statistics 2024/2025 _________________________________

#### 6.1.GEOSTATISTICS ####
#__________________________#

## Clear the workspace
rm(list=ls())

## Load spatial packages

library(sp)           ## Data management
library(lattice)      ## Data management
library(gstat)        ## Geostatistics (essential package)
library(geoR)         ## Geostatistics (just for some plots)
#                        (for people with Mac: you need to have XQuartz https://www.xquartz.org/)

## Set working directory 

## Functions for graphics 
v.f <- function(x, ...){100-cov.spatial(x, ...)}
v.f.est<-function(x,C0, ...){C0-cov.spatial(x, ...)}



####   EXPLORATORY ANALYSIS & VARIOGRAM ESTIMATION  ####

## Load meuse data set:
## The meuse is a classical geostatistical data set used frequently
## to demonstrate various geostatistical analysis steps.
## The point data set consists of 155 samples of top soil heavy metal
## concentrations (ppm), along with a number of soil and landscape variables.
## The samples were collected in a flood plain of the river Meuse,
## near the village Stein (The Netherlands).

data(meuse) # available in gstat
?meuse

# before starting working, be sure you have transformed latitude & longitude
# in UTM coordinates (planar coordinates, in meters)

## Define the sample coordinates
coordinates(meuse) <- c('x','y')
meuse # meuse is now a spatial dataframe

# bubble plot(obj,zcol,...)
# key.space=location of the key
bubble(meuse, 'zinc', do.log=TRUE, key.space='bottom')
# dimension of of the bubbles is proportional to the value of variable 'zinc'
# --> we notice that the higher values of zinc are close to river

spplot(meuse, "zinc", do.log = T, colorkey = TRUE)

dev.off()


# just for a better visualization of the river Meuse:
# we load river meuse data
data(meuse.riv)
meuse.lst <- list(Polygons(list(Polygon(meuse.riv)), "meuse.riv"))
meuse.sr <- SpatialPolygons(meuse.lst)
# grid for prediction
data(meuse.grid)
is(meuse.grid)
coordinates(meuse.grid) <- c('x','y')
meuse.grid <- as(meuse.grid, 'SpatialPixelsDataFrame')
# plot all together
image(meuse.grid, col = "lightgrey")
plot(meuse.sr, col = "grey", add = TRUE)
plot(meuse, add = TRUE)
title('meuse river geostatistical data')
# sampling locations are crosses
# we visualize the river
# domain D is grey area (subset of R^2 because we used UTM coordinates)

dev.off()


##### Exploratory Analysis #####
# we select a variable, e.g. zinc
# [as homework, you could replicate the analysis with other variables]

# histogram of zinc variable 
hist(meuse$zinc, breaks=16, col="grey", main='Histogram of Zn', prob = TRUE, xlab = 'Zn')
# highly skewed, transform to the log
hist(log(meuse$zinc), breaks=16, col="grey", main='Histogram of log(Zn)', prob = TRUE, xlab = 'log(Zn)')

# scatterplot of log(zinc) with respect to distance from the river 
xyplot(log(zinc) ~ sqrt(dist), as.data.frame(meuse))
# Negative Correlation: lower distance from the river => higher level of zinc
# it supports what we observed in the bubble plot
dev.off()



##### Estimating Spatial Correlation #####
#           Variogram Analysis     

# sample (empirical) variogram (binned estimator)
# the function variogram, which takes a formula as its first argument: 
#.      the notation "~ 1" stands for a single constant predictor
#.      (hp: spatially constant mean - stationary model)
# log(zinc)~1 means that we assume a constant trend for the variable log(zinc).

svgm <- variogram(log(zinc) ~ 1, meuse) # with all the defaults of the function

plot(svgm, main = 'Sample Variogram', pch=19)
# Observations:
# 1) is there an asymptote? if yes, it is compatible with stationarity;
# 2) behaviour near zero? Linear model could fit well; also a quadratic might work.

# 15 lags by default
# default decisions:
# cutoff, lag width, direction dependence

# cutoff distance: maximum distance up to which point pairs are considered
#                  (default: bbox diagonal / 3) [empirical rule]
#                  where bbox stands for bounding box
# NB: cutoff distance should be much shorter than the maximum distance among points

# lag width: width of distance intervals over which point pairs are averaged
#            in bins (default: cutoff distance / 15)

plot(variogram(log(zinc) ~ 1, meuse, cutoff = 1000, width = 1000/15),pch=19)  # smaller cutoff
plot(variogram(log(zinc) ~ 1, meuse, cutoff = 1000, width = 1000/5),pch=19)

# intervals can have different widths: to fix varying widths use the argument boundaries
c(0,200, seq(400,1500,100))
plot(variogram(log(zinc) ~ 1, meuse, boundaries = c(0,200, seq(400,1500,100))),pch=19)
# change it depending on the disposition of your data in space


# the following
plot(variogram(log(zinc) ~ 1, meuse), pch=19)
# automatically decides to ignore direction: point pairs are merged on the
# basis of distance to compute the empirical variogram

# we can take into account different directions through 'alpha':
plot(variogram(log(zinc) ~ 1, meuse, alpha = c(0, 45, 90, 135)), pch=19)
# point pairs whose separation vector has a given direction 
# are used in each panel (not too many directions otherwise noise will increase)
# Note: Zonal Anisotropy 
# look at direction 45 --> difference in the asymptote --> sill --> less variability
# it is the direction parallel to the river --> we should include the river in the model

# for the moment we proceed assuming isotropy, 
# but, in general, the options we have are the following:
# --> model anisotropy (anisotropic variogram model; it can be done by the package)
# --> correct anisotropy by including a non-stationarity,
#     accounting for anisotropy in a covariate



##### Variogram modeling #####

# we go on assuming isotropicity

# list of parametric isotropic variogram models
vgm()

# basic variograms available
show.vgms() 

# in gstat, valid variogram models are constructed by using one or
# combination of two or more basic variogram models
# first argument of the function 'vgm' is partial sill,
# then the desired model, then range, and finally nugget: 
# vgm(sill, model, range, nugget)

# some examples...
vgm(1, "Sph", 300)
vgm(1, "Sph", 300, 0.5)

# one can also add two or more models
v1 <- vgm(1, "Sph", 300, 0.5)
v2 <- vgm(0.8, "Sph", 800, add.to = v1) # by specifying "add.to = .. " 
v2

# this is only measurement error
vgm(0.5, "Nug", 0)


## weighted least squares fitting a variogram model to the sample variogram
## STEPS:
## 1) choose a suitable model
## 2) choose suitable initial values for partial sill, range & nugget
## 3) fit the model using one of the possible fitting criteria

v <- variogram(log(zinc) ~ 1, meuse)
plot(v,pch=19)

# Linear behavior near the origin, growth not very fast 
# Recall: both spherical and exponential model have a Linear behavior near the
#         origin but Exponential Model has a Faster Growth than the Spherical One
# => we fit a spherical model

# try reasonable initial values (remember: vgm(sill, model, range, nugget))
fit.variogram(v, vgm(0.65, "Sph", 800, 0.1))
fit.variogram(v, vgm(1, "Sph", 800, 1)) # --> sill and nugget are largely overestimated;
                                        # the most difficult parameter to fit is the range

# try unreasonable initial values
fit.variogram(v, vgm(1, "Sph", 10, 1))
# due to high non-linearity in the minimization problem,
# starting from unreasonable initial values might cause fail to converge

# plot of the final fit
v <- variogram(log(zinc) ~ 1, meuse)
v.fit <- fit.variogram(v, vgm(1, "Sph", 800, 1))
plot(v, v.fit, pch = 19) # always do this!

# fitting method: non linear regression with minimization of weighted
# sum of squares error. Final value of the minimum
attr(v.fit, 'SSErr')
# how can we choose weights? argument fit.method in fit.variogram:
# fit.method = 1 : w = N_j
# fit.method = 2 : w = N_j/gamma(h_j)^2
# fit.method = 6 : w = 1
# fit.method = 7 : w = N_j/h_j^2 (default method)

# one can also keep one of the parameters fixed, and fit only the others.
# This is common for the nugget parameter, which may be hard to infer from data
# when sample locations are regularly spread. Information may be derived from
# measurement error characteristics for a specific device.

# ex: fix the nugget variance to the value 0.06
fit.variogram(v, vgm(1, "Sph", 800, 0.06), fit.sills = c(FALSE, TRUE))
# the range parameters can be fixed using argument fit.ranges

## maximum likelihood fitting of variogram models
## - does not need the sample variogram
## - can be performed through restricted maximum likelihood
fit.variogram.reml(log(zinc) ~ 1, meuse, model=vgm(0.6, "Sph", 800, 0.06))
# compare it with
v.fit


## modeling anisotropy*
v.dir <- variogram(log(zinc)~1, meuse, alpha=(0:3)*45) 
v.anis <- vgm(.6, "Sph", 1600, .05, anis=c(45, 0.3))

print(plot(v.dir, v.anis, pch=19))
# --> better to correct anisotropy by including a non-stationarity,
#     accounting for anisotropy in a covariate



####          SPATIAL PREDICTION & KRIGING          ####

##### Stationary Univariate Spatial Prediction (Ordinary Kriging) #####

## Prediction at a single new location 
s0.new = data.frame(x=179180, y=330100) # UTM coordinates 
coordinates(s0.new) = c('x','y')        # set the coordinates of the dataframe

# plot all together
image(meuse.grid, col = "lightgrey")
plot(meuse.sr, col = "grey", add = TRUE)
plot(meuse, add = TRUE)
plot(s0.new, add = TRUE, col='red', lwd = 2)
title('meuse river geostatistical data')

# Create a gstat object setting a spherical (residual) variogram
# gstat(g.obj, id, formula, data, model, set,...)
g.tr <- gstat(formula = log(zinc) ~ 1, data = meuse, model = v.fit) # v.fit is the variogram model we use
g.tr

## ordinary kriging
# Make the ordinary kriging prediction with the function: 
# predict(obj, grid, BLUE=FALSE)
# this gives the prediction of Y(s_0):
predict(g.tr, s0.new)
# var1.pred is the predicted log(zinc) 
# var1.var variance of prediction error (ordinary kriging variance)
# variance > 0 (as expected)

# Estimate the mean: use the argument 'BLUE'
predict(g.tr, s0.new, BLUE = TRUE)
# this gives the estimate (best linear unbiased estimator) of the mean
# (trend component) under gls

## consider a location where I OBSERVE data
# this gives the prediction of Y(s_0)
# in the first location (zero variance!)
meuse[1,]
log(1022)
predict(g.tr, meuse[1,]) 
# --> prediction is exactly that value
# --> zero prediction variance: I know that point

# this gives the estimate of the mean
# (drift component) under gls
predict(g.tr, meuse[1,], BLUE = TRUE) # same mean everywhere (I am under stationarity, same as before)

# prediction over the entire grid
lz.ok <- predict(g.tr, meuse.grid, BLUE = FALSE)


spplot(lz.ok[1], main='Ordinary Kriging - Prediction')
spplot(lz.ok[2], main='Ordinary Kriging - Variance')


##### Non-stationary Univariate Spatial Prediction (Universal Kriging) #####

# the hypothesis of spatially constant mean may be too restrictive!
# we now use as covariate the square root of the distance from the river Meuse

# to fit the variogram on the residuals, one should take into account 
# the spatial dependence while estimating the trend component by using GLS

# Create a gstat object setting a spherical (residual) variogram
# gstat(g.obj, id, formula, data, model, set,...)
meuse.gstat <- gstat(id = 'zinc', formula = log(zinc) ~ sqrt(dist), # non stationary formula (intercept is automatically inside)
                     data = meuse, nmax = 50, model = v.fit, set = list(gls=1))
# nmax = 50 is the maximum number of iterations (50 is actually much larger than needed)
# model = v.fit is the initial model used for iterations (here we just give the one estimated under stationarity)
# set = list(gls=1) means that I want the estimation to be done with GLS
meuse.gstat

# Estimate the variogram from GLS residuals:
?variogram.gstat
v.gls <- variogram(meuse.gstat) # difference wrt before (before we were just giving a formula, now we give the non-stationary object)
plot(v.gls) 
# before it was 0.6, now 0.25

v.gls.fit <- fit.variogram(v.gls, vgm(0.25, "Sph", 800, 0.8)) # vgm(sill, model, range, nugget)
v.gls.fit
plot(v.gls, v.gls.fit, pch = 19)

# Update gstat object with variogram model
meuse.gstat <- gstat(id = 'zinc', formula = log(zinc) ~ sqrt(dist),
                     data = meuse, nmax = 50, model = v.gls.fit, set = list(gls=1))
meuse.gstat


## universal kriging:
s0.new
## I have to define also the covariate in s_0
s0.vec <- as.vector(slot(s0.new,'coords'))
# distance to the river: calculate the distance between s0 and each point of
# the river, then select the minimum
s0.dist <- min(rowSums(scale(meuse.riv,s0.vec)^2)) 
s0.new <- as.data.frame(c(s0.new,s0.dist))
names(s0.new) <- c('x','y','dist')
coordinates(s0.new) <- c('x','y')
s0.new <- as(s0.new, 'SpatialPointsDataFrame')
s0.new

# Function "predict" uses the residual variogram stored in the gstat
# object to make the prediction
predict(meuse.gstat, s0.new) # you need to give object and new location
# variance > 0 (as expected) --> large underestimation of uncertainty (--> sigma is actually unknown)

# this gives the estimate of x(s_0)'*beta
# (trend component) under gls
predict(meuse.gstat, s0.new, BLUE = TRUE)

# prediction over the entire grid
lz.uk <- predict(meuse.gstat, meuse.grid, BLUE=FALSE)

# estimate of the mean over the entire grid
lz.uk.BLUE <- predict(meuse.gstat, meuse.grid, BLUE=TRUE)

spplot(lz.ok[,1], main = 'Ordinary Kriging') # pattern is given by the data

x11()
spplot(lz.uk[,1], main = 'Universal Kriging') # much stronger patterns following the river

x11()
spplot(lz.uk.BLUE[,1], main = 'Universal Kriging - drift')

x11()
spplot(lz.uk[,2], main = 'Universal Kriging - Variance')


# Is the drift important to explain the variability of the response variable z_s?
# z_s = m_s (drift) + delta_s (residuals)
# Let's compare the variogram of the data and of the residuals:
plot(v$dist,v$gamma,xlab='distance',ylab='semivariance',pch=19,col='skyblue1',ylim=c(0,0.8))
curve(v.f.est(x, C0=v.fit[2,2]+v.fit[1,2], cov.pars=rbind(c(v.fit[2,2], v.fit[2,3]),c(v.fit[1,2], v.fit[1,3])), 
              cov.model = c("spherical","pure.nugget")), from = 0.0001, to = 1600,
              xlab = "distance", ylab = expression(gamma(h)),
      main = "Variogram model",add=TRUE,col='skyblue1',lwd=2, ylim=c(0,110)) # variance of z estimated from observations (0.6)

points(v.gls$dist,v.gls$gamma,xlab='distance',ylab='semivariance',pch=19,col='steelblue',ylim=c(0,0.8))
curve(v.f.est(x, C0=v.gls.fit[2,2]+v.gls.fit[1,2], 
              cov.pars=rbind(c(v.gls.fit[2,2], v.gls.fit[2,3]),c(v.gls.fit[1,2], v.gls.fit[1,3])), 
              cov.model = c("spherical","pure.nugget")), from = 0.0001, to = 1600,
      xlab = "distance", ylab = expression(gamma(h)),
      main = "Variogram model",add=TRUE,col='steelblue',lwd=2, ylim=c(0,110)) # variance of the residuals (0.2) [iterative algorithm]

# variance explained by the drift is what is in between the two curves;
# the larger the difference between the sills of stationary and non stationary model, 
# the more the regressors are important in explaining variability
# if there is not much difference among the two, we can go for a stationary model (simpler)
graphics.off()



# xx <- variogram(log(zinc) ~ sqrt(dist), data = meuse)
# fit.xx <- fit.variogram(xx, vgm(1, "Sph", 800, 1))
# yy <- gstat(formula = log(zinc) ~ sqrt(dist), data = meuse, model = fit.xx)
# predict(yy, s0.new, BLUE = TRUE)



# Further references on this example
# - documentation of the package
#   https://www.gstat.org/gstat.pdf
# - a brief tutorial
#   https://cran.r-project.org/web/packages/gstat/vignettes/gstat.pdf
# - further references with extra plots, functions and materials
#   http://statweb.lsu.edu/faculty/li/IIT/spatial.html
#   https://rstudio-pubs-static.s3.amazonaws.com/558732_e09c44c118d94efb951ab71b26abf80b.html
#   http://rstudio-pubs-static.s3.amazonaws.com/10213_8c02d102993942a88574e44abdf3a235.html
