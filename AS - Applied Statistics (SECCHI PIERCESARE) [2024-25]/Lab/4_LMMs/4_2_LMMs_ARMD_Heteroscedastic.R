#_________________ Applied Statistics 2024/2025 _________________________________

#### 4.3.Linear Mixed-effects models: ARMD Trial - Heteroscedastic Residuals ####
#_______________________________________________________________________________#

# Topics:
#   LINEAR MIXED MODELS WITH HETEROSCEDASTIC RESIDUALS (VarPower())
#   3. Linear Models with random intercept 
#   4. Linear Models with random intercept + slope 
#      4.A general structure of D
#      4.B diagonal D 
#--> Interpretation of random effects
#--> Diagnostic
#--> Models comparison 
#--> Tests for the random components



#### LINEAR MIXED MODELS WITH HETEROSCEDASTIC RESIDUALS ####


# fixed-effects formula
lm2.form <- formula(visual ~ visual0 + time + treat.f + treat.f:time ) 

# LMM with homoscedastic residuals and random intercept
fm16.1 <- lme(lm2.form, random = ~1|subject, data = armd)

#________________________________________________________________________________________
#### Model 3. Random intercept, heteroscedastic residuals (varPower of time) (M16.3) ####


# update fm16.1 including heteroscedastic residuals
fm16.3 <- update(fm16.1,
                 weights = varPower(form = ~ time), 
                 data = armd) # using REML

summary(fm16.3)
# sigma = 3.607 (residual standard deviation)
# delta = 0.3144
# sigma * sqrt(d11) = 7.706 (standard deviation of the random intercepts)

VarCorr(fm16.3)  

## var-cov matrix of the errors (i.e. of Y, conditional to the random effects), 
# that are independent but heteroscedastic 
# sigma^2 * R_i
fm16.3ccov = getVarCov(fm16.3, type = "conditional",  individual = "2")
fm16.3ccov
# example: sigma^2 * 4^(2*delta) = 3.60672^2 * 4^(2*0.3144) = 31.103

plot(as.matrix(fm16.3ccov[[1]]), 
     col=colorRampPalette(c("blue", "white", "red")),
     main = expression(paste('Conditional estimated Var-Cov matrix of ', Y[2])))


## var-cov matrix of Y_i
# sigma^2 * V_i
fm16.3cov = getVarCov(fm16.3, type = "marginal", individual = "2")
fm16.3cov # (90.479 = 31.103 + 59.37555; 121.440 = 62.062 + 59.37555; ...)
# The estimated variance of random intercepts sigma^2*d11 = 59.376. 
# --> It is smaller than the value of 80.608, obtained for model 16.1
# This is expected, because, by allowing for heteroscedastic
# residual random errors, a larger part of the total variability is explained by the
# residual variances.

# Corr(sigma^2 * V_i)
cov2cor(fm16.3cov[[1]])                    
# The corresponding estimated marginal correlation matrix indicates a decreasing
# correlation between visual acuity measurements made at more distant timepoints.

plot(as.matrix(fm16.3cov[[1]]), 
     col=colorRampPalette(c("blue", "white", "red")),
     main = expression(paste('Marginal estimated Var-Cov matrix of ', Y[2])))

# NB. Direct comparison of the marginal variance covariance
# matrices for models M12.3 and M16.2 is not appropriate. 
# --> marginal variance-covariance matrix of model M16.2, 
# is much more structured than that of model M12.3. 
# On the other hand, they both allow for marginal correlation coefficients, 
# which depend on the time “distances”, or “positions”, 
# of visual-acuity measurements.




## ANALYSIS OF RESIDUALS

# Default residual plot of conditional Pearson residuals
plot(fm16.3)
# the plot is not very informative,
# because it pools all the residuals together, despite the fact that residuals 
# obtained from the same individual are potentially correlated. 
# However, it can serve for detecting, e.g., outliers.
# (at the bottom and the top of the central part of the scatterplot)


# Boxplots of Pearson residuals per time and treatment
bwplot(resid(fm16.3, type = "p") ~ time.f | treat.f, 
       data = armd)
# Despite standardization, the variability of the residuals seems to vary.
# The plot reveals also a number of outliers (in all treatment groups and at all timepoints)


# Normal Q-Q plots of Pearson residuals 
qqnorm(fm16.3, ~resid(.) | time.f) 


## ANALYSIS OF RANDOM EFFECTS
# Normal Q-Q plots of predicted random effects
qqnorm(fm16.3, ~ranef(.))  
# NB observed distribution of bi^ does not necessary reflect the true distribution of bi 
# --> interpret with caution


## Computing predictions comparing population average predictions with patient-specific predictions
aug.Pred <- augPred(fm16.3,
                    primary = ~time, # Primary covariate
                    level = 0:1,     # fixed/marginal (0) and subj.-spec.(1)
                    length.out = 2)  # evaluated in two time instants (4 e 52 wks)

plot(aug.Pred, layout = c(4, 4, 1))
# The predicted population means decrease linearly in time.
# According to the assumed structure of the model, the population
# means are shifted for individual patients by subject-specific random intercepts.
# Slopes of the individual profiles are the same for all subjects.

# For example, for the subjects 4 and 15, the
# predicted individual patterns suggest a decrease of visual acuity over time, while the
# observed values actually increase over time.
# A possible way to improve the individual predictions is to allow not only for
# patient-specific random intercepts, but also for patient-specific random slopes.


#________________________________________________________________________________________________________________
#### Model 4.A. random intercept + slope (correlated), heteroscedastic residuals (varPower of time) (M16.4A) ####

fm16.4A <- update(fm16.3,
                 random = ~1 + time | subject,
                 data = armd)
summary(fm16.4A)

getVarCov(fm16.4A, individual = "2")  # sigma^2 * D_i (i=2)

intervals(fm16.4A, which = "var-cov")  # Estimates of sigma^2 * D, delta e sigma

# The results show a low estimated value
# of the correlation coefficient for the random effects b0i and b2i, equal to 0.138.
# The confidence interval for the correlation coefficient suggests that, in fact, the two
# random effects can be uncorrelated. 
# --> we consider a simplified form of the D matrix.


#______________________________________________________________________________________________________#
## Model 4.B. random intercept + slope independent, heteroscedastic residuals (varPower of time) (M16.4B)
fm16.4B <- update(fm16.4A,
                 random = list(subject = pdDiag(~time)), # Diagonal D
                 data = armd) 
summary(fm16.4B) ## results suggest to remove the Treat and Time interaction

intervals(fm16.4B)
# delta = 0.11 --> marginal variance function is quadratic over time
# the var(Yit) function increases with time (since d11, d22 and sigma^2 are positive)


anova(fm16.4B, fm16.4A)  # H0: d_12=0
# We test if d_12 = 0 --> not statistically significant at 5%, we can simplify the D structure in diagonal


## ANALYSIS OF RESIDUALS
qqnorm(fm16.4B, ~resid(.) | time.f) 

## ANALYSIS OF RANDOM EFFECTS
qqnorm(fm16.4B, ~ranef(.)) 
# to be interpreted with caution since it might not reflect the real unknown distribution

bwplot(resid(fm16.4B, type = "p") ~ time.f | treat.f, 
       panel = panel.bwplot, # User-defined panel (not shown)
       data = armd)

## We make predictions comparing population average predictions with patient specific predictions
aug.Pred <- augPred(fm16.4B,
                    primary = ~time, # Primary covariate
                    level = 0:1,     # Marginal(0) and subj.-spec.(1)
                    length.out = 2)  # Evaluated in two time instants (4 e 52 wks)
plot(aug.Pred, layout = c(4, 4, 1), columns = 2) 
# slopes vary. -->  the predicted individual profiles follow more closely 
# the observed values and capture, e.g., increasing trends in time,
# --> better fit than before


## let's compare the 4 fitted models 
AIC(fm16.1,fm16.3,fm16.4A, fm16.4B)
anova(fm16.1,fm16.3,fm16.4A, fm16.4B)



