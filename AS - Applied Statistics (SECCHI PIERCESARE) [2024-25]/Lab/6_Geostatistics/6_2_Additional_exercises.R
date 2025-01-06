#_________________ Applied Statistics 2024/2025 _________________________________

#### 6.2.Additional exercises ####
#________________________________#


##### exercise 1) #####

# One of the most relevant consequences of the eruption of volcan 
# Eyjafjoll (in Iceland), in 2010, is the contamination by fluoride. 
# The latter is due to the deposit on the ground of the ash released
# in the atmosphere during the eruption.
# The file "fluoruro.txt" reports the coordinates of 50 measurement sites
# s_i, i=1,...,50, the corresponding concentrations of fluoride (ppm) F(s_i)
# and the distance D.s_i of each site s_i from the crater of the volcano.
# Denoting by delta a zero-mean, second-order stationary and isotropic random
# field:
# a) Estimate two empirical variograms, assuming the following models:
#    F(s_i)=beta_0+delta(s_i) and 
#    F(s_i)=beta_0+beta_1*D.s_i+delta(s_i). 
#    Choose the most appropriate model for the observations.
# b) Fit to the empirical variogram at point (a), a Gaussian model
#    without nugget, via weighted least squares. Use as initial parameters: 
#    sill=100, range=0.02. Report the estimates of sill and range.
# c) Fit to the empirical variogram chosen at point (a), a spherical model
#    without nugget, via weighted least squares. Report the estimates of sill 
#    and range.
# d) Compare the variograms estimated at points (b) and (c), with the empirical
#    variogram at point (a). Given that the ash deposition is known to be
#    a very regular phenomenon, which variogram model is the most appropriate?
# e) Based on model (d), estimate the concentration of fluoride due to the eruption
#    in the city of Raufarhofn (s0 = (0.3; 0.24), D.s0 = 0.1970) 
# f) Based on model (d), estimate the concentration of fluoride at the same
#    location, due to a possible new eruption of equivalent intensity, independent 
#    of that of 2010.

#  SOLUTION
#_____________

# Data import

library(sp)          ## Data management
library(lattice)      ## Data management
library(gstat)        ## Geostatistics (essential package)

data=read.table('fluoruro.txt')
names(data)[3]='f'
attach(data)
coordinates(data)=c('X','Y')
head(data) # --> correctly imported

# a) Estimate two empirical variograms, assuming the following models:
#    F(s_i)=beta_0+delta(s_i) and 
#    F(s_i)=beta_0+beta_1*D.s_i+delta(s_i). 
#    Choose the most appropriate model for the observations.

v=variogram(f ~ 1, data=data)   # stationary model: do we like it?
plot(v,pch=19)

v.t=variogram(f ~ D, data=data) # non stationary model --> better --> we clearly see a range & sill
plot(v.t,pch=19)

# b) Fit to the empirical variogram at point (a), a Gaussian model
#    without nugget, via weighted least squares. Use as initial parameters: 
#    sill=100, range=0.02. Report the estimates of sill and range.
v.fit2 <- fit.variogram(v.t, vgm(100, "Gau", 0.02))
plot(v.t, v.fit2, pch = 3)
v.fit2

# c) Fit to the empirical variogram chosen at point (a), a Spherical model
#    without nugget, via weighted least squares. Report the estimates of sill 
#    and range.
v.fit1 <- fit.variogram(v.t, vgm(100, "Sph", 0.08))
plot(v.t, v.fit1, pch = 3)
v.fit1

# d) Compare the variograms estimated at points (b) and (c), with the empirical
#    variogram at point (a). Given that the ash deposition is known to be
#    a very regular phenomenon, which variogram model is the most appropriate?

# We choose v.fit2 --> we now that ash deposition is very regular
# and fit2 generates smoother prediction because quadratic in the origin


# e) Based on model (d), estimate the concentration of fluoride due to the eruption
#    in the city of Raufarhofn (s0 = (0.3; 0.24), D.s0 = 0.1970) 

g.t <- gstat(formula = f ~ D, data = data, model = v.fit2)

D.s0=0.1970
s0=as.data.frame(matrix(c(0.3,0.24,D.s0),1,3))
names(s0)=c('X','Y','D')
coordinates(s0)=c('X','Y')

predict(g.t, s0, BLUE = FALSE)

# f) Based on model (d), estimate the concentration of fluoride at the same
#    location, due to a possible new eruption of equivalent intensity, independent 
#    of that of 2010.

predict(g.t, s0, BLUE = TRUE) # best guess of the mean of the process
# new eruption is totally independent with the eruption of 2010
# I am back to the case of independence
# I just know something on the mean of the new process

detach(data)



#________________________________________________________________________________
#

##### exercise 2) #####

# The file radioville.txt reports the information on 158 control units
# in the area around the nuclear power plant of Radioville.
# At each site, available data consist of: radioactivity levels [Bq],
# longitude [°N], latitude [°W] and type of soil [urban/vegetation].
# Denoting by s_i the i-th site, by R the radioactivity level,
# by eps a weakly stationary random field and by D a dummy urban/vegetation:

# a) estimate the parameters of the linear model 
#    R(s_i) = beta_0 + beta_1 D(s_i) + eps(s_i) 
#    assuming for eps a spherical variogram without nugget, estimated
#    via weighted least squares;
# b) estimate the parameters of the linear model 
#    R(s_i) = beta_0 + beta_1 D(s_i) + eps(s_i) 
#    assuming for eps a spherical variogram with nugget, estimated
#    via weighted least squares;
# c) choose the best variogram model by comparing the fitted model 
#    with the corresponding empirical variograms (report qualitative
#    plots and the estimated variogram parameters)
# d) on the basis of model (c), predict the radioactivity level at the
#    parking lot of the shopping centre of Radioville (lon = 78.59, 
#    lat = 35.34), and in the park of Radioville (lon = 77.6, 
#    lat = 34.99);
# e) estimate variance of prediction error at the same locations
#    as at point d).

#  SOLUTION
#_____________

data <- read.table('radioville.txt',header=TRUE)
attach(data)

# create dummy: 0 = urban, 1 = vegetation
DUMMY <- rep(0,length(D))
DUMMY[which(D=='V')] <- 1
data <- data.frame(cbind(Bq,Long,Lat,DUMMY))
names(data) <- c('Bq','Long','Lat','D')
coordinates(data) <- c('Long','Lat')
data

## point a)
## fitting a variogram without nugget

v <- variogram(Bq ~ D, data = data) # empirical semivariogram
plot(v)
v.fit1 <- fit.variogram(v, vgm(1, "Sph", 0.5))
plot(v, v.fit1, pch = 3)
v.fit1

# coefficient of the linear model: 
# it sufficies to estimate the drift at two locations where we have observations,
# with D=U and D=V
# data[1,] = urbane        (D=0)
# data[6,] = vegetation    (D=1)
g.tr <- gstat(formula = Bq ~ D, data = data, model = v.fit1)
predict(g.tr, data[1,], BLUE = TRUE)
predict(g.tr, data[6,], BLUE = TRUE)

## point b)
## fitting a variogram with nugget

v <- variogram(Bq ~ D, data = data)
plot(v)
v.fit2 <- fit.variogram(v, vgm(0.6, "Sph", 0.5, 0.1))
plot(v, v.fit2, pch = 3)
v.fit2

# it sufficies to estimate the drift at two locations where we have observations,
# with D=U and D=V
# data[1,] = urbane
# data[6,] = vegetation
g.tr <- gstat(formula = Bq ~ D, data = data, model = v.fit2)
predict(g.tr, data[1,], BLUE = TRUE)
predict(g.tr, data[6,], BLUE = TRUE)


## point d)
## predict at 2 new locations: we use model 1 (without nugget)
g.tr <- gstat(formula = Bq ~ D, data = data, model = v.fit1)

# urbane : 78.59,35.34
s0.new <- as.data.frame(matrix(c(78.59,35.34,0),1,3))
names(s0.new) <- c('lon','lat','D')
s0.new
coordinates(s0.new) <- c('lon','lat')
predict(g.tr, s0.new)

# vegetation : 77.69,34.99
s0.new <- as.data.frame(matrix(c(77.69,34.99,1),1,3))
names(s0.new) <- c('lon','lat','D')
s0.new
coordinates(s0.new) <- c('lon','lat')
predict(g.tr, s0.new)


