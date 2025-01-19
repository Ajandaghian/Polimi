#__________ Applied Statistics 2024/2025 _____________

#### 3.3.LMs with correlated errors: Hawaii Birds ####
#____________________________________________________#

# Source: Mixed Effects Models and Extensions in Ecology with R (2009) Zuur, Ieno, Walker, Saveliev and Smith. Springer

# Annual abundances of three bird species measured at three islands in Hawaii from 1956 to 2003.
setwd("/Users/amirh_jandaghian/Documents/Polimi/Y 02/Polimi - Git/AS - Applied Statistics (SECCHI PIERCESARE) [2024-25]/Lab/3_LMs_fixed_effects_correlated_errors")
# We focus on Moorhen.Kauai:
Hawaii = read.table('Hawaii.txt', header=T)
Hawaii$Birds <- sqrt(Hawaii$Moorhen.Kauai)
plot(Hawaii$Year, Hawaii$Moorhen.Kauai, xlab = "Year", ylab = "Moorhen abundance on Kauai")
# general increase since the mid 1970s.

library(nlme)
lm0 = gls(Moorhen.Kauai ~ Rainfall + Year, na.action = na.omit, data = Hawaii)
summary(lm0)
plot(lm0) # Pearson residuals are plotted

# problem with the homoscedasticity --> heterogeneous variance

# we try to stabilize it by applying a square root transformation
lm1 = gls(sqrt(Moorhen.Kauai) ~ Rainfall + Year, na.action = na.omit, data = Hawaii)
plot(Hawaii$Year, sqrt(Hawaii$Moorhen.Kauai), xlab = "Year", ylab = "Moorhen abundance on Kauai")
plot(lm1)

# or with the logarithm
lm2 = gls(log(Moorhen.Kauai) ~ Rainfall + Year, na.action = na.omit, data = Hawaii)
plot(Hawaii$Year, log(Hawaii$Moorhen.Kauai), xlab = "Year", ylab = "Moorhen abundance on Kauai")
plot(lm2)

# take home message: sometimes it is enough to transform the features;
# in this example we keep it simple and we go on in this way, instead of introducing a variance function

Hawaii$Birds <- sqrt(Hawaii$Moorhen.Kauai)
M0 <- gls(Birds ~ Rainfall + Year, na.action = na.omit, data = Hawaii)
summary(M0)
# The summary table shows that the effect of Rainfall is not significant, 
# but there is a significant increase in birds over time. 
# --> we cannot trust these p-values as we may be violating the independence assumption. 


# Birds_s = beta0 + beta1 * Rainfall_s + beta2 * Year_s + eps_s
# eps_s ∼ N(0,sigma^2)
# cov(eps_s, eps_t) = {sigma^2 if s=t
#                     {0       otherwise

# --> we model the auto-correlation between residuals of different time points by introducing a function h():
# cor(eps_s, eps_t) = {1                    if s=t
#                     {h(eps_s, eps_t, rho) otherwise

# We assume stationarity, 
# e.g. correlation between the residuals eps_s and eps_t only depends on their time difference s – t. 
# Hence, the correlation between eps_s and eps_t is assumed to be
# the same as that between eps_s+1 and eps_t+1, between eps_s+2 and eps_t+2, etc.

# --> we need to model h()


# To detect patterns:
# --> standardized residuals
plot(M0)

# --> auto-correlation function (ACF)
# The value of the ACF at different time lags gives an indication whether there is any auto-correlation in the data
E <- residuals(M0, type = "normalized")
I1 <- !is.na(Hawaii$Birds)
Efull <- vector(length = length(Hawaii$Birds))
Efull <- NA
Efull[I1] <- E   # we insert the missing values that were removed by gls with na.omit
acf(Efull, na.action = na.pass, main = "Auto-correlation plot for residuals")

# --> clear violation of the independence assumption; 
# --> various time lags have a significant correlation!
# --> the ACF plot has a general pattern of decreasing values for the first 5 years


# Under the independence assumption, C is a identity matrix;
# We introduce now:

# 1) Compound Symmetry Structure:
#    It assumes that whatever the distance in time between two observations,
#    their residual correlation is the same
# cor(eps_s, eps_t) = {1       if s=t
#                     {rho     otherwise

M1 <- gls(Birds ~ Rainfall + Year,
          na.action = na.omit, data = Hawaii ,
          correlation = corCompSymm(form =~Year))
summary(M1)
anova(M0, M1) # --> no improvements


# 2) AR-1 auto-correlation.
# cor(eps_s, eps_t) = {1              if s=t
#                     {rho^(|t-s|)    otherwise

M2 <- gls(Birds ~ Rainfall + Year,
          na.action = na.omit, data = Hawaii,
          correlation = corAR1(form =~ Year))
summary(M2)
anova(M0,M2)  # --> improvement

# Occasionally, a negative rho can be found
# Plausible explanation is the abundances go from high values in one year 
# to low values in the next year.

# 3) ARMA Error Structures
# ARMA model has two parameters defining its order: 
# the number of auto-regressive parameters (p) 
# and the number of moving average parameters (q)
# if p=1 and q=0 we have AR-1 model
# in general, better not to select a model with p,q > 2 or 3 (convergence problems)
cs1 <- corARMA(c(0.2), p = 1, q = 0) # we provide arbitrary starting points
cs2 <- corARMA(c(0.3, -0.3), p = 2, q = 0)
M3arma1 <-gls(Birds ~ Rainfall + Year,
              na.action = na.omit,
              correlation = cs1, data = Hawaii)
M3arma2 <- gls(Birds ~ Rainfall + Year,
               na.action = na.omit,
               correlation = cs2, data = Hawaii)
AIC(M0, M2, M3arma1, M3arma2)

summary(M3arma2)
# Phi1 close to 1 may indicate a more serious problem of the residuals 
# being non-stationary (non-constant mean or variance).
# Note that the auto-correlation function in the plot becomes positive again for
# larger time lags, suggesting that an error structure that allows for a sinusoidal
# pattern may be more appropriate.



# NB: "there is not much to be gained from finding the perfect correlation structure 
# compared to finding one that is adequate".
# [Schabenberger and Pierce (2002), Diggle et al. (2002), andVerbeke and Molenberghs (2000)]



