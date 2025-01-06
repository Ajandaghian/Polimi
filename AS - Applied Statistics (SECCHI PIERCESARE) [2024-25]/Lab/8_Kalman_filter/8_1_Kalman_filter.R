#_________________ Applied Statistics 2024/2025 _________________________________

#### 8.1. Kalman filter  ####
#___________________________#

# install.packages("KFAS")
library(KFAS)

# Sources: 
# https://www.jstatsoft.org/article/view/v078i10
# https://cran.r-project.org/web/packages/KFAS/vignettes/KFAS.pdf

# Our time series consists of yearly alcohol-related deaths per 100,000 persons 
# in Finland for the years 1969–2007 in the age group of 40–49 years

data("alcohol", package = "KFAS")
deaths <- window(alcohol[, 2], end = 2007)
population <- window(alcohol[, 6], end = 2007)

ts.plot(deaths / population, 
        ylab = "Alcohol-related deaths in Finland per 100,000 persons", 
        xlab = "Year", col = c(1))


#### Example of Gaussian state space model (Section 2.2) ####

Zt <- matrix(c(1, 0), 1, 2)
Ht <- matrix(NA) # unknown variance parameter to be estimated with fitSSM
Tt <- matrix(c(1, 0, 1, 1), 2, 2)
Rt <- matrix(c(1, 0), 2, 1)
Qt <- matrix(NA) # unknown variance parameter to be estimated with fitSSM
a1 <- matrix(c(1, 0), 2, 1)
P1 <- matrix(0, 2, 2)
P1inf <- diag(2)


model_gaussian <- SSModel(deaths / population ~ -1 + 
                            SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf),
                          H = Ht)

# SSModel function is the formula which defines 
# - the observations (left side of tilde operator ~) 
# - and the structure of the state equation (right side of tilde)

# deaths / population is a univariate time series,
# state equation is defined using the system matrices with auxiliary function SSMcostum
# the intercept term is omitted with -1 in order to keep the model identifiable

##### gaussian_fit #####

# The unkown variance parameters can be estimated using the function fitSSM

fit_gaussian <- fitSSM(model_gaussian, inits = c(0, 0), method = "BFGS")

# maximum likelihood estimates
## \sigma^2_{\epsilon}
fit_gaussian$model$H
## \sigma^2_{\eta}
fit_gaussian$model$Q

out_gaussian <- KFS(fit_gaussian$model)
# these are the final estimates
out_gaussian 

## constant slope term
out_gaussian$a[40, 2]
## with standard error
sqrt(out_gaussian$P[2, 2, 40])

##### gaussian_plot #####
ts.plot(cbind(deaths / population, 
              out_gaussian$a[-(length(deaths) + 1), 1], out_gaussian$alphahat[, 1]), 
        ylab = "Alcohol-related deaths in Finland per 100,000 persons", 
        xlab = "Year", col = c(1, 2, 4))

# observations with one-step-ahead predictions (red) and smoothed (blue)
# estimates of the random walk process µ_t
# At the time t the Kalman filter computes the one-step-ahead prediction error vt = yt − µt,
# and uses this and the previous prediction to correct the prediction for the next time point
# this is most easily seen at the beginning of the series,
# where our predictions seem to be lagging the observations by one time step.
# The smoothing algorithm takes account of both the past and the future values at each
# time point, thus producing more smoothed estimates of the latent process


#### Example of non-Gaussian state space model (Section 3.2) ####

# Alcohol-related deaths can also be modeled naturally as a Poisson process.
# Now our observations yt are the actual counts of alcohol-related deaths in year t, 
# whereas the varying population size is taken into account by the exposure term ut

# we now need to define the distribution of
# the observations using the argument distribution (which defaults to "gaussian"). 
# We also define the exposure term via the argument u 
# (for non-Gaussian models the H is omitted and vice versa), 
# and use default values for a1 and P1 in the SSMcustom.

model_poisson <- SSModel(deaths ~ -1 + 
                           SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, P1inf = P1inf), 
                         distribution = "poisson", 
                         u = population)

##### nongaussian_fit #####
fit_poisson <- fitSSM(model_poisson, inits = -5, method = "BFGS")

out_poisson <- KFS(fit_poisson$model, smoothing = "state")
out_poisson

## \sigma^2_{\eta} (the only unknown parameter)
fit_poisson$model$Q
# but the actual values between the Gaussian and Poisson models are not directly comparable
# as the interpretation of µt differs between models.

## slope term of the Poisson model
out_poisson$alphahat[1, 2]
## with standard error
out_poisson$V[2, 2, 1]
## % yearly increase in deaths
signif(100 * (exp(out_poisson$alphahat[1, 2]) - 1), 2)


##### nongaussian_plot #####
x11()
ts.plot(deaths / population, out_gaussian$alphahat[, 1], 
        exp(out_poisson$alphahat[, 1]), xlab = "Year",
        ylab = "Alcohol-related deaths in Finland per 100,000 persons", 
        col = c(1,4,2))

# estimates of the intensity (deaths per 100,000 persons) modeled
# as Gaussian process (blue), and as a Poisson process (red).



#### Example of ARIMA model (Section 6.2) ####

# we again model the alcohol-related deaths but now use the ARIMA(0, 1, 1)
# model with drift:
  
##### arima_time_series #####

drift <- 1:length(deaths)

model_arima <- SSModel(deaths / population ~ drift + 
                         SSMarima(ma = 0, d = 1, Q = 1))
# The auxiliary function SSMarima defines the ARIMA model

update_model <- function(pars, model) {
  tmp <- SSMarima(ma = pars[1], d = 1, Q = pars[2])
  model["R", states = "arima"] <- tmp$R
  model["Q", states = "arima"] <- tmp$Q
  model["P1", states = "arima"] <- tmp$P1
  model
}

# In this case we need to supply the model updating function for fitSSM which updates our
# model definition based on the current values of the parameters we are estimating. 
# Instead of manually altering the corresponding elements of the model, 
# update_model uses the SSMarima function for computation of relevant 
# system matrices R, Q and P1

fit_arima <- fitSSM(model_arima, inits = c(0, 1), updatefn = update_model, 
                    method = "L-BFGS-B", lower = c(-1, 0), upper = c(1, 100))
fit_arima$optim.out$par

## \theta_1, \sigma
c(fit_arima$model["R", 4],
  fit_arima$model["Q"])


