#__________________________ Applied Statistics 2024/2025 _________________________

#### 2.3.Squid Dataset: more on LMs with heteroscedastic and independent errors ####
#_________________________________________________________________________________#

# Data published by Smith et al. (2005), 
# who looked at seasonal patterns in reproductive and somatic tissues 
# in the squid 'Loligo forbesi'.
# --> we'll focus on dorsal mantle length (in mm) and testis weight from 768 male squid.

# The aim is to model the testis weight as a function of the dorsal mantle length (DML)
# and the month recorded, 
# for determining the degree to which sexual maturation is size-related and seasonal.

# Source: Mixed Effects Models and Extensions in Ecology with R (2009) 
# Zuur, Ieno, Walker, Saveliev and Smith. Springer

Squid = read.table('Squid.csv', header = TRUE, dec=".", sep=',')

str(Squid)
Squid$fMONTH = as.factor(Squid$MONTH)

#_________________________________________________________________________________
# let's fit our first "classical" model

# Testisweight_im = beta_0 + beta_1*DML_i + beta_2m + beta_3m * DML_i + eps_im,  
# for i=1,...,768, m=1,...,12
# with beta_2m = beta_3m = 0 if m=1

# eps_im ~ N(0,sigma^2)
M1 <- lm(Testisweight ~ DML * fMONTH, data = Squid)
summary(M1)
anova(M1)
# all regression parameters are significantly different from 0 at the 5% level

# Automatic diagnostic plots (graphically)
x11()
par(mfrow = c(2,2))
plot(M1)
shapiro.test(M1$residuals)

# we focus on this plot:
plot(M1, which = c(1), col = 1, add.smooth = FALSE, caption = "")
# Clear violation of homoscedasticity: 
# the larger the length (DML), the larger the variability;

# We try to detect why:
# we look at MONTH
colori = rainbow(length(unique(Squid$fMONTH)))
x11()
boxplot(M1$residuals ~ Squid$fMONTH, col=colori,
        xlab='Months', ylab='Residuals', 
        main ='Distribution of residuals across months')  
# or more easily: plot(Squid$fMONTH, resid(M1), xlab = "Month", ylab = "Residuals")

# we look at DML
x11()
plot(Squid$DML, resid(M1, type = "pearson"), xlab = "DML", ylab = "Residuals")

# we cannot trust these results as we are clearly violating the homogeneity assumption:
# the assumption that var(eps_i) = sigma^2 is wrong
# --> we might assume that the residuals have variance var(eps_i) increasing with DML_i

#_______________________________________________________________________________

# R function	   Description
# VarFixed	     Fixed variance
# VarIdent	     Different variances per stratum
# VarPower	     Power of the variance covariate
# VarExp	       Exponential of the variance covariate
# VarConstPower	 Constant plus power of the variance covariate
# VarComb        A combination of variance functions

#_______________________________________________________________________________

# Testisweight_im = beta_0 + beta_1*DML_i + beta_2m + beta_3m * DML_i + eps_im,  
# for i=1,...,768, m=1,...,12
# with beta_2m = beta_3m = 0 if m=1

# eps_im ~ N(0, sigma^2 * DML_i) for i=1,...,768
# Such a variance structure allows for larger residual spread if DML increases.
# No extra parameters involved.

library(nlme)
M.lm <- gls(Testisweight~DML*fMONTH, data=Squid)
# NB. unweighted least square is a classical linear model

# Fixed variance along DML (--> variance that is proportional to DML).
vf1Fixed <- varFixed(~DML)
M.gls1 <- gls(Testisweight~DML*fMONTH, weights=vf1Fixed, data=Squid)

anova(M.lm,M.gls1)
# models are not nested; 
# --> no LR test statistic is given but AIC favours model M.gls1!
# we can also use: 
AIC(M.lm, M.gls1)

# but...
plot(M.gls1, which = c(1), col = 1, add.smooth = FALSE, caption = "")

#_______________________________________________________________________________

# Now allow for different variance for each month
# eps_im ~ N(0, sigma_m^2) for i=1,...,768 and m=1,...,12
vf2 <- varIdent(form= ~ 1|fMONTH) # different variances per stratum for month
M.gls2 <- gls(Testisweight ~ DML*fMONTH, weights=vf2, data =Squid)
summary(M.gls2)

anova(M.lm,M.gls2)
# we are testing whether
# H0: sigma_1^2 = sigma_2^2 = sigma_3^2 = ... = sigma_12^2
# vs
# H1: they are different

# we have evidence to reject the null hypothesis that all variances are the same: 
# model with different variances per month is better.
# Notice the df: due to the variance structure, we now have to estimate 11 more parameters. 

M.gls2$modelStruct$varStruct  # one multiplication factor is set to 1 (in this case, month2)
summary(M.gls2)$sigma
# In month 9, the variance is (2.99 * 1.27)^2, in month 12 it is (1.27 * 1.27)^2
# Note that months 9, 10 and 3 have the highest ratio, 
# indicating that in these months there is more residual variation.


# So, which option is better: different spread per MONTH or different spread along DML?
boxplot(M.lm$residuals ~ Squid$fMONTH, col=colori,
        xlab='Months', ylab='Residuals', 
        main ='Distribution of residuals across months')  
plot(M.lm, which = c(1), col = colori, add.smooth = FALSE, caption = "")
# the smaller fitted values are from months with less spread
# and the larger fitted values are from months with higher spread

# Let’s look at the residuals across months to further investigate.
library(lattice)
x11()
coplot(resid(M.lm)~DML|fMONTH,data=Squid)
# The lower left panel corresponds to month 1, the lower right to month 4, and
# the upper right to month 12;
# The residual variation differs per month, but in some
# months (e.g. 3, 9, and 10) the residual spread also increases for larger DML values. 
# So, both are influential: residual spread is influenced by both month and length!


# Before discussing how to combine both types of variation 
# (variation linked with DML and variation linked with Month), 
# we introduce a few more variance structures.

# The easiest approach for choosing the best variance structure 
# is to apply the various available structures in R and compare them 
# using the AIC or to use biological knowledge combined with some 
# integrative graphs like the coplot.

# Some of the variance functions are nested, and a LR test can be applied to 
# judge which one performes better for the data.


#_______________________________________________________________________________

##### varPower #####

# eps_im ~ N(0, sigma^2 * |DML_im|^(2*delta)) for i=1,...,768 and m=1,...,12
# delta is unknown and needs to be estimated;

# If delta=0   we obtain the "classical" linear regression model,
#              and the two models are nested 
#              --> LR test can be applied to judge which one is better. 
# if delta=0.5 and variance covariate with positive values:
#              we get the same variance structure as varFixed(~DML). 
# If variance covariate has values equal to 0, the variance of the residuals is 0 as well. 
# This causes problems in the numerical estimation process and varPower should not be used. 
# [For the squid data, all DML values are larger than 0 (DML is length)]

vf3 <- varPower(form =~ DML)
M.gls3 <- gls(Testisweight ~ DML * fMONTH, weights = vf3, data=Squid)
summary(M.gls3)                # --> delta=1.759
AIC(M.lm,M.gls1,M.gls2,M.gls3) # lowest value so far


# we can also model an increase in spread for larger DML values, but only in certain months!
vf4 <- varPower(form=~ DML | fMONTH)
M.gls4 <- gls(Testisweight ~ DML * fMONTH, data = Squid, weights = vf4)
AIC(M.lm,M.gls1,M.gls2,M.gls3,M.gls4) # lowest value so far
summary(M.gls4)                       # we can get the values of the deltas (we have 12)
# little variation, but they are multiplied by two before being used to take the power


# It is also possible to set the delta_j for some months equal to 
# an a priori chosen value and keep it fixed. 
# --> e.g. if you know or want to test whether the spread along DML in some months is constant
vf4_f <- varPower(form=~ DML | fMONTH, fixed=c('4'=1.75, '11'=1.75))
M.gls4_f <- gls(Testisweight ~ DML * fMONTH, data = Squid, weights = vf4_f)
AIC(M.lm,M.gls1,M.gls2,M.gls3,M.gls4,M.gls4_f) # lowest value so far
summary(M.gls4_f)                              # we can get the values of the deltas (we have 12)


#_______________________________________________________________________________
##### varExp() #####

# eps_im ~ N(0, sigma^2 * exp(2 * delta * |DML_im|) for i=1,...,768 and m=1,...,12
# delta is unknown and needs to be estimated;
# --> if delta = 0, gives the "classical model"
# --> there are no restrictions on delta or DML
# --> this structure allows also a decrease of spread for DML values if delta is negative.
vf5 <- varExp(form =~ DML)
M.gls5 <- gls(Testisweight ~ DML * fMONTH,
              weights = vf5, data = Squid)
summary(M.gls5)
AIC(M.lm,M.gls1,M.gls2,M.gls3,M.gls4,M.gls5)  # slightly higher than M.gls3

# As before, we can allow for different delta per month ( varExp(form =∼ DML | fMONTH) )
# Again, it is possible to fix some of the delta

#_______________________________________________________________________________

##### varConstPower() #####
# constant plus power of the variance covariate function

# eps_im ~ N(0, sigma^2 * (delta_1 + |DML_im|^delta_2)^2 ) for i=1,...,768 and m=1,...,12
# If delta_1 = 1 and delta_2 = 0, we are back to the "classical" linear regression model
# If not, then the variance is proportional to a constant plus the power of the variance covariate DML.

# This variance structure should work better than the varExp if the variance
# covariate has values close to zero.
vf6 <- varConstPower(form =~ DML)
M.gls6 <- gls(Testisweight ~ DML * fMONTH, weights = vf6, data = Squid)
summary(M.gls6)
AIC(M.lm,M.gls1,M.gls2,M.gls3,M.gls4,M.gls5,M.gls6)


# Again, we can allow for different delta1s and delta2s per stratum of a nominal variable 
# (e.g. MONTH).
# eps_im ~ N(0, sigma^2 * (delta_1m + |DML_im|^delta_2m)^2 ) for i=1,...,768 and m=1,...,12
vf7 <- varConstPower(form =~ DML | fMONTH)
M.gls7 <- gls(Testisweight ~ DML * fMONTH, weights = vf7, data = Squid)
AIC(M.lm,M.gls1,M.gls2,M.gls3,M.gls4,M.gls5,M.gls6,M.gls7)
# Again, it is possible to set the delta1s and delta2s to a preset value for
# particular months and keep it fixed during the estimation process.

#_______________________________________________________________________________

##### varComb() #####

# Combination of variance structures using the varComb function
# --> we can allow for both an increase in residual spread for larger DML values 
#     as well as a different spread per month.
# For example, we want to combine varIdent and varExp.
# eps_im ~ N(0, sigma_m^2 * exp(2*delta*DML_im) for i=1,...,768 and m=1,...,12
# sigma has an index m running from 1 to 12, allowing for different spreads per month.
vf8 <- varComb(varIdent(form =~ 1 | fMONTH), varExp(form =~ DML) )
M.gls8 <- gls(Testisweight ~ DML * fMONTH, weights = vf8, data = Squid)

M.gls8$modelStruct$varStruct 

anova(M.lm,M.gls1,M.gls2,M.gls3,M.gls4,M.gls5,M.gls6,M.gls7,M.gls8)
# The model allowing for an increase in spread for larger DML values 
# (which is allowed to differ per month), M.gls4, has the lowest AIC and is therefore
# selected as the optimal model.

# seems that varPower(form=~ DML | fMONTH) is the best model
anova(M.lm, M.gls4)

#_______________________________________________________________________________

# NB.
# If the variance covariate has large values (e.g. larger than 100), 
# numerical instabilities may occur; exp(100) is rather large! 
# In such cases, it is better to rescale the
# variance covariate before using it in any of the variance structures. For example, we
# could have used DML/max(DML) or express it in meters instead of millimetres in
# the variance functions. The unscaled DML can still be used in the fixed part of the
# model.

# NB. Nested structures
# - "Classical" linear model is not nested within varFixed structure
#   because we have eps_i~N(0,sigma^2) VS eps_i~N(0, sigma^2 * DML_i)
#   We cannot obtain the variance structure on the left from the right one, 
#   unless DML is equal to 1 for all observations. 
# - "Classical" linear model is nested within varPower structure
#   because we have eps_i~N(0,sigma^2) VS eps_i~N(0, sigma^2 * |DML_ij|^(2*delta_j))
#   By setting all delta_j equal to zero in the right variance structure
# - varIdent is nested in the varPower structure

# NB. Which variance structure should you choose?
# --> If the variance covariate is a categorical variable, use varIdent.
# --> varFixed, varPower, varExp, and varConstPower  allow for an increase (or decrease) in
#     residual variation for the Testisweight data along a continuous variance covariate
#     like DML (an explanatory variable in this case)
# --> varFixed is rather limited, as it assumes that the variance of the residuals 
#     is linearly related to a variance covariate. This causes problems if the variance
#     covariate takes non-positive values or where the linear relationship requirements
#     between variation and the variance covariate is too stringent.
# --> it may be better to use the varPower, varExp, or varConstPower functions, 
#     which allow for more flexibility.
# --> The varPower should not be used if the variance covariate takes the value of zero.
# --> VarConstPower structure should work better than varExp if the variance covariate 
#     has values close to zero.


#_______________________________________________________________________________

# Let’s graphically confirm.
# First look at ordinary residuals
E1 <- resid(M.gls4)
x11()
coplot(E1 ~ DML | fMONTH,ylab="Ordinary residuals", data = Squid)

# Note that these residuals still show heterogeneity, 
# but this is now allowed (because the residual variation differs
# depending on the chosen variance structure and values of the variance covariate).
# Hence, these residuals are less useful for the model validation process.

# Let's look now at the normalized residuals
# eps_im = (Testisweight_im - Fitted values_im) / sqrt(sigma^2 * |DML_im|^(2*delta_m))

E2 <- resid(M.gls4, type = "normalized")
x11()
coplot(E2 ~ DML | fMONTH, data = Squid, ylab = "Normalised residuals") # homoscedasticity

anova(M.gls4)




