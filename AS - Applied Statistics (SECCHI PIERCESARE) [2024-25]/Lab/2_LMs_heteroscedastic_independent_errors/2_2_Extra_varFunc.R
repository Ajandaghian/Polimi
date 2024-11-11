#______________________________ Applied Statistics 2024/2025 ________________________________

#### 2.2.Extra about varFunc class (from Chapter 8) ####
#_______________________________________________________


?Initialize.varFunc
# parameters of Initialize:
# Initialize( object of varFunc class, data )

# Example with 
(val <- c("12wks" = 0.5, "24wks" = 2))  # delta1 = 1, delta2 = 0.5, delta3 = 2     
(fix <- c("52wks" = 3))                 # delta4 = 3 (fixed)

# the initial value of SD at 12 weeks is specified as a half of that at 4weeks. 
# the value of SD at 24 weeks and 52 weeks is taken as twice and three times as high, 
# respectively, as the value at 4weeks. 
# The use of the argument fixed=fix implies that the value of the coefficient corresponding to
# the variance for week 52 will not change during any optimization steps in modeling routines.

frm  <- formula(~1|time.f)              # time.f is a stratifying factor
(vf0 <- varIdent(value = val,           # Var. function object defined... 
                 fixed = fix,
                 form  = frm)) 
(vf0i <- Initialize(vf0, armd))         # ... and initialized


# (Chapter 8.3)
# parameters of coef:
# coef( object of class varClass, 
#       parametrization applied to the coefficients:unconstrained?, 
#       which coefficients
#     )
# where:
# - unconstrained (logical): the coefficients (parameters) can be presented on a constrained 
#   or unconstrained scale (cfr 7.4.3). (default: unconstrained = TRUE)
# - allCoef (logical): indicates whether all coefficients or only those, 
#   which were not designated to be fixed in numerical optimization routines, are to be returned.
#   (default: allCoef = FALSE)


coef(vf0i, unconstrained = FALSE, allCoef = TRUE) # All delta coefs
coef(vf0i, unconstrained = FALSE, allCoef = FALSE)# Varying only

coef(vf0i, unconstrained = TRUE, allCoef = TRUE)  # All delta* coefs
coef(vf0i, unconstrained = TRUE, allCoef = FALSE) # Varying (default)

# NB. for objects of class varIdent(), coefficients on the unconstrained scale 
# are obtained by taking the natural logarithm of the corresponding constrained parameters

coef(vf0i) <- c(-0.6, 0.7)                        # New coefs assigned   
coef(vf0i, allCoef = TRUE)                        # All coefs printed


# Useful syntax: extracting information from a varFunc class object.
summary(vf0i)               # Summary
formula(vf0i)               # Variance function formula
getCovariate(vf0i)          # Variance covariate
getGroupsFormula(vf0i)      # Formula for variance strata
stratum <-getGroups(vf0i)   # Stratification variable
length(stratum)             # Length of stratum indicator
unique(stratum)             # Unique strata
stratum[1:6]                # First six observations
varWeights(vf0i)[3:6]       # Variance weights 1/^lambda_i
logLik(vf0i)                # Contribution to the log-likelihood
