#_____________ Applied Statistics 2024/2025 _________________________

#### 3.2.Linear Models with heteroscedastic and correlated errors ####
#____________________________________________________________________#

# Topic:
# Linear Model with Fixed Effects and Correlated Errors (theory: Chapter 10)
# Chapters 11 and 12 (R) + ARMD dataset

library(nlme)
library(lattice)


#### The gls() function ####

#gls.fit = gls(
#              model =
#              data =
#              subset =
#              na.action =
#              method =
#              weights = [used to model the Variance Function, saw in lab2]
#              correlation = corStruct(form=formula)    [default = NULL]
#          )


# generic functions to extract information: print(), summary() and results().
#_Component to be extracted_______________Syntax________________________________________
# gls()-call                              (cl <- getCall(gls.fit))
# correlation = argument                  cl$correlation
# 95% CI for rho                          intervals(gls.fit, which ="var-cov")$corStruct
# ^R_i matrices                           getVarCov(gls.fit)
# Normalized residuals                    resid(gls.fit, type = "normalized")
# Var-cov structure                       mSt <- gls.fit$modelStruct
# Correlation structure (CorSt)           cSt <- mSt$corStruct

# Summary                                 summary(cSt)
# CorSt formula                           formula(cSt)
# CorSt covariate                         getCovariate(cSt)
# Contribution to log-likelihood          logLik(cSt)
# ^C_i matrices                           corMatrix(cSt)  



#corStruct class
?corClasses

# (Chapter 11.2)
# corStruct class specifies correlation structures for the model-fitting function gls()
# with the help of the "correlation" argument.

# Generally, objects of the class corStruct have 3 arguments:
# - 'value' : correlation parameter;
# - 'form' : one sided formula (position variable and optionally a grouping factor);
#            observations in different groups are assumed to be uncorrelated;
#            default: form = ~1
# - 'fixed' : TRUE or FALSE; if TRUE, values of correlation parameters are fixed in numerical opt


#_______________________________________________________________________________
# In addition (but we will see these things later in spatial statistics)
# For Spatial Correlation functions (e.g. corExp, corGaus,..) there is also:
# - 'nugget' : FALSE or TRUE; if FALSE: 'value' assumes only one value (the correlation parameter, aka range>0)
#                             if TRUE: 'value' = c(range, nugget); we recall that range>0 and 0<nugget<1
#                                      default is value=numeric(0): i.e nugget=0.1 and range=0.9*(min of pairwise distances);
# - 'metric' = "euclidean" (default) or "maximum" or "manhattan"

# NB. for spatial corr functions, form = ~s1+...+sP|g, where
#     - s1...sP are spatial position variables (can be unidimensional and we go back to the longitudinal data);
#     - g is an optional grouping factor. When it is present, we are assuming correlation only
#       between some level of grouping factor: different levels of observations are uncorrelated.











