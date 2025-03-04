#_________________ Applied Statistics 2024/2025 _______________________________

#### 4.1.Linear Mixed-effects models: ARMD Trial - Homoscedastic Residuals ####
#_____________________________________________________________________________#

library(nlmeU)
library(corrplot)
library(nlme)
library(lattice)
library(plot.matrix)
library(lme4)
library(insight)

rm(list=ls())
graphics.off()

# Topics:
#   LINEAR MIXED MODELS WITH HOMOSCEDASTIC RESIDUALS
#   1. Linear Models with random intercept 
#   2. Linear Models with random intercept + slope 
#      2.A general structure of D
#      2.B diagonal D 
#--> Interpretation of random effects and PVRE
#--> Diagnostic
#--> Models comparison 
#--> Tests for the random components


data(armd, package = "nlmeU") # Age-Related Macular Degeneration
rlims <- c(-4.5, 4.5)
xlims <- c(0, 90)
xlimc <- c("4", "12", "24", "52wks")


## Linear Mixed-Effects models (LMM)
## In the LMM approach, the hierarchical structure of the data is directly addressed, 
## with random effects that describe the contribution of the variability at different levels 
## of the hierarchy to the total variability of the observations.

## Two main R packages: 
## 1. 'lme4' with lmer() --> it does not handle heteroscedastic residuals
## 2. 'nlme' with lme()  --> it handles heteroscedastic residuals

## --> We will mainly use lme() function in nlme package for LMM models;
## --> We will shift to lmer() function in lme4 package for LMM models with homoscedastic residuals (this script)
##     only for some visualizations



# LINEAR MIXED MODELS WITH HOMOSCEDASTIC RESIDUALS 

#___________________________________________________________________
#### Model 1. Random intercept, homoscedastic residuals (M16.1) ####

# We now treat time as a numeric variable (as assessed in the last lab)
lm2.form <- formula(visual ~ visual0 + time + treat.f + treat.f:time)

(fm16.1 <- lme(lm2.form, random = ~1|subject, data = armd))   
# By default, lme() assumes independent residual errors with a constant variance, sigma^2.
# The default REML estimation is used --> to change it to the ML estimation, put method="ML" 

summary(fm16.1)
# sigma*sqrt(d11) = 8.98 (standard deviation of the random intercepts)
# sigma = 8.63           (residual standard deviation)


# print out the estimated fixed-effects table
printCoefmat(summary(fm16.1)$tTable, has.Pvalue = TRUE, P.values = TRUE) 
# Remember: t-test statistics for the fixed-effects coefficients (marginal-approach)

# to get the confidence intervals of all the estimated parameters
intervals(fm16.1)


# Var-Cov matrix of fixed-effects
vcovb <- vcov(fm16.1) 
vcovb
# and Correlation of fixed effects
corb <- cov2cor(vcovb) 
nms <- abbreviate(names(fixef(fm16.1)), 5)
rownames(corb) <- nms
corb


# Var-Cov matrix of random-effects and residuals (what we observe in the summary)
print(vc <- VarCorr(fm16.1), comp = c("Variance", "Std.Dev."))
VarCorr(fm16.1)
var_eps = as.numeric(vc[2,1])
var_eps

sd_eps <- summary(fm16.1)$sigma
sd_eps

var_b = as.numeric(vc[1,1])
var_b


# Let's compute the conditional and marginal var-cov matrix of Y
#---------------------------------------------------------------#
# --> Conditional variance-covariance matrix of Y given b (diagonal matrix) 
#     sigma^2 * Ri for the second subject (type='conditional')

getVarCov(fm16.1,                     
          type = "conditional",       # sigma^2 * R_i = 74.434 * I4
          individual = "2")   
# Conditioned to the random effects b_i --> var-cov of the errors are independent and homoscedastic

# we extract sigma^2 * Ri for patients i=2,3,4
sR = getVarCov(fm16.1, type = "conditional", individual = 2:4)
# and we plot them
plot(as.matrix(bdiag(sR$`2`, sR$`3`, sR$`4`)), 
     col=colorRampPalette(c("blue", "white", "red")),
     main = 'Conditional estimated Var-Cov matrix of Y given b')



# --> the marginal variance-covariance matrix of Y (block-diagonal matrix)
#     sigma^2 * V_i for the second subject (type='marginal')
(sVi = getVarCov(fm16.1,                       
          type = "marginal",      # sigma^2 * V_i: sigma^2*d11 extra-diagonal and sigma^2+d11 on main diagonal
          individual = "2"))
(cov2cor(sVi[[1]]))               # Corr(sigma^2 * V_i) 

# we extract sigma^2 * V_i for patients i=2,3,4
sV <- getVarCov(fm16.1, type = "marginal", individual = 2:4) 
# visualization 
plot(as.matrix(bdiag(sV$`2`, sV$`3`, sV$`4`)), #-> V is a block-diagional matrix, the marginal var-cov matrix
     col=colorRampPalette(c("blue", "white", "red")),
     main = 'Marginal estimated Var-Cov matrix of Y')


# PVRE 
#------#
# i.e. the Percentage of Variance explained by the Random Effect (PVRE).
# This is also called the intraclass correlation (ICC),
# because it is also an estimate of the within group correlation.

PVRE <- var_b/(var_b+var_eps)
PVRE # it is high!



# Visualization of confidence intervals
#--------------------------------------#
intervals(fm16.1)

# fixed effects
# we plot the "classical" linear model for visualizing the fixed effects (without intercept)
fm16.1GLS <- gls(lm2.form, data = armd) 
library(sjPlot)
plot_model(fm16.1GLS)
plot_model(fm16.1)

# random effects
## visualization of the random intercepts with their 95% confidence intervals
# Random effects: b_0i for i=1,...,234
re = ranef(fm16.1)
dat = data.frame(x= row.names(re),y=re[,attr(re,'effectName')])
# The dotplot shows the point and interval estimates for the random effects
# ordered
dotplot(reorder(x,y)~y,data=dat)
# not ordered
plot(ranef(fm16.1))

#_____________________________________________________________________________
# nicer visualization of the random effects with lmer() function
# formulation with lmer() for including a random intercept at subject level
fm16.1mer <- lmer(visual ~ visual0 + time * treat.f + (1|subject), data = armd)
# we can highligh which are significantly different from the mean (0)
# visualization of the random intercepts with their 95% confidence intervals
dotplot(ranef(fm16.1mer)) 

#install.packages('TMB', type = 'source')
plot_model(fm16.1mer, type='re') #--> positive (blu) and negative (red) effect
#_____________________________________________________________________________




#Test for (variance of) random intercept
#--------------------------------------#
# H0 is "variance of a random effect is 0 
#       in an LMM with a known correlation structure of the tested random effect and
#       independent and identically distributed random errors"
library(RLRsim)    
exactRLRT(fm16.1)        # M16.1 (alternative model) 
# The function simulates values of the REML-based LRtest statistic.

# The p-value of the REML-based LR test,
# estimated from 10,000 simulations (the default), 
# --> pval clearly indicates the importance of including the random intercepts into the model, 
#     needed to adjust for the correlation between visual acuity measurements.

# Because we test a random effect in model M16.1,
# which contains only a single random effect (only the intercept), 
# we use the abbreviated form of the function call, with model M16.1 as the only argument. 


# Diagnostic plots 
#----------------#
# 1) Assessing Assumption on the within-group errors
plot(fm16.1)  # Pearson and raw residuals are the same now 
# (no scale is applied since we are dealing with homogeneous variance)

qqnorm(resid(fm16.1)) # normality of the residuals
qqline(resid(fm16.1), col='red', lwd=2)

# 2) Assessing Assumption on the Random Effects
qqnorm(fm16.1, ~ranef(.), main='Normal Q-Q Plot - Random Effects on Intercept')



#_______________________________________________________________________________
#### Model 2: random intercept + slope and homoscedastic residuals  (M16.2) ####

##### Model 2.A: general D (M16.2A) #####
#---------------------------------------#

fm16.2A <- lme(lm2.form, random = ~1 + time | subject, data = armd)

summary(fm16.2A)
intervals(fm16.2A)

# variance covariance matrix of the fixed effects
vcovb <- vcov(fm16.2A) 
vcovb
# and correlation
corb <- cov2cor(vcovb) 
nms <- abbreviate(names(fixef(fm16.2A)), 5)
rownames(corb) <- nms
corb

# Var-Cov matrix of random-effects and errors
print(vc <- VarCorr(fm16.2A), comp = c("Variance", "Std.Dev."))
sigma <- summary(fm16.2A)$sigma


# Let's compute the conditional and marginal var-cov matrix of Y
#---------------------------------------------------------------#
# the conditional variance-covariance matrix of Y (diagonal matrix)
getVarCov(fm16.2A,                     
          type = "conditional",       # sigma^2 * R_i
          individual = "2")   

# we extract sigma^2 * R_i for patients i=2,3,4
sR = getVarCov(fm16.2A, type = "conditional", individual = 2:4)
# and we plot them
plot(as.matrix(bdiag(sR$`2`, sR$`3`, sR$`4`)), 
     col=colorRampPalette(c("blue", "white", "red")),
     main = 'Conditional estimated Var-Cov matrix of Y given b')

# the marginal estimated variance-covariance matrix of Y (block-diagonal matrix)
(sVi = getVarCov(fm16.2A,                       
                type = "marginal",      # sigma^2 * V_i
                individual = "2"))
(cov2cor(sVi[[1]]))                     # Corr(sigma^2 * V_i) 

# we extract sigma^2 * V_i for patients i=2,3,4
sV <- getVarCov(fm16.2A, type = "marginal", individual = 2:4)
# and we plot them
plot(as.matrix(bdiag(sV$`2`, sV$`3`, sV$`4`)), 
     col=colorRampPalette(c("blue", "white", "red")),
     main = 'Marginal estimated Var-Cov matrix of Y')


# PVRE 
#-----#
# In this case the variance the mean random effect variance of the model is given by
# var_b = Var(b0,b1) = sigma2_b0 +  sigma2_b1 * mean(w^2) + 2 * Cov(b0,b1) * mean(w) 
# remember that Cov(b0,b1) = Cor(b0,b1) * sd(b0) * sd(b1)
# See equation (10) in Johnson (2014), Methods in Ecology and Evolution, 5(9), 944-946.

var_eps <- as.numeric(vc[3,1])
var_eps
var_b <- as.numeric(vc[1,1]) + as.numeric(vc[2,1])*mean(armd$time^2) + 
  + 2*as.numeric(vc[2,3])*as.numeric(vc[1,2])*as.numeric(vc[2,2])*mean(armd$time) 
var_b


PVRE <- var_b/(var_b+var_eps)
PVRE # 72% is very high!



# Visualization 
#-------------#
# visualization of the random intercepts & slopes with their 95% confidence intervals
# Random effects: b_0i, b_1i for i=1,...,234
plot(ranef(fm16.2A))


# as before, we can fit better plots through lmer()
fm16.2Amer <- lmer(visual ~ visual0 + time * treat.f + (1+time|subject), data = armd, 
                   control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
dotplot(ranef(fm16.2Amer))

# for a better visualization
plot_model(fm16.2Amer, type='re') #--> positive (blu) and negative (red) effect



# Comparing models 
#-----------------#
# The anova function, when given two or more arguments representing fitted models,
# produces likelihood ratio tests comparing the models.
anova(fm16.1, fm16.2A)
# The p-value for the test is essentially zero -> we prefer fm16.2A



# Diagnostic plots
#-----------------#
# 1) Assessing Assumption on the within-group errors
plot(fm16.2A)

qqnorm(resid(fm16.2A))
qqline(resid(fm16.2A), col='red', lwd=2)

# 2) Assessing Assumption on the Random Effects
qqnorm(fm16.2A, ~ranef(.), main='Normal Q-Q Plot - Random Effects on Intercept & Slope')



intervals(fm16.2A, which = "var-cov")
# We observed that the 95% CI for correlation between sigma^2*d_11 and sigma^2*d_22 id contains 0 --> 
# --> they can be uncorrelated
# --> we fit a new model with a diagonal D matrix 



#_______________________________________________________________________________
##### Model 2.B: diagonal D (M16.2B) ######
#-------------------------------------------------------------------------------#

fm16.2B <- lme(lm2.form, random = list(subject = pdDiag(~time)), data = armd) 
# Diagonal D (diagonal positive-definite matrix)

intervals(fm16.2B)                       # 95% CI for betas, sigma  

summary(fm16.2B)

vcovb <- vcov(fm16.2B) 
vcovb
corb <- cov2cor(vcovb) 
nms <- abbreviate(names(fixef(fm16.2B)), 5)
rownames(corb) <- nms
corb

# Var-Cov matrix of random-effects and residuals
print(vc <- VarCorr(fm16.2B), comp = c("Variance", "Std.Dev."))


sigma <- summary(fm16.2B)$sigma
sigma


# Let's compute the conditional and marginal var-cov matrix of Y
#---------------------------------------------------------------#
# the conditional variance-covariance matrix of Y (diagonal matrix)
getVarCov(fm16.2B,                     
          type = "conditional",       # sigma^2 * R_i
          individual = "2")   
# Conditioned to the random effects b_i, we observe the var-cov of the errors
# that are independent and homoscedastic

# we extract sigma^2 * R_i for patients i=2,3,4
sR = getVarCov(fm16.2B, type = "conditional", individual = 2:4)
# and we plot them
plot(as.matrix(bdiag(sR$`2`, sR$`3`, sR$`4`)), 
     col=colorRampPalette(c("blue", "white", "red")),
     main = 'Conditional estimated Var-Cov matrix of Y given b')

# the marginal variance-covariance matrix of Y (block-diagonal matrix)
(sVi = getVarCov(fm16.2B,                       
                type = "marginal",      # sigma^2 * V_i
                individual = "2"))
(cov2cor(sVi[[1]]))                     # Corr(sigma^2 * V_i) 

# we extract sigma^2 * V_i for patients i=2,3,4
sV <- getVarCov(fm16.2B, type = "marginal", individual = 2:4)
# and we plot it
plot(as.matrix(bdiag(sV$`2`, sV$`3`, sV$`4`)), 
     col=colorRampPalette(c("blue", "white", "red")),
     main = 'Marginal estimated Var-Cov matrix of Y')


# PVRE
#-----#
# In this case the variance of random effects represents the mean random 
# effect variance of the model and is given by
# var_b = Var(b0,b1) = sigma2_b0 + 0 + sigma2_b1*mean(z^2)
# See equation (10) in Johnson (2014), Methods in Ecology and Evolution, 5(9), 944-946.
vc
var_eps <- as.numeric(vc[3,1])
var_eps
var_b <- as.numeric(vc[1,1]) + mean(armd$time^2)*as.numeric(vc[2,1]) # 54.07117 + 0.07935904*mean(armd$time^2) 
var_b

PVRE <- var_b/(var_b+var_eps)
PVRE # 72% is very high!


# Visualization
#--------------#
# visualization of the random intercepts with their 95% confidence intervals
# Random effects: b_0i, b_1i for i=1,...,234
plot(ranef(fm16.2B))

# for a better visualization
fm16.2Bmer <- lmer(visual ~ visual0 + time * treat.f + (1|subject) + (0 + time|subject),
                   data = armd, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
dotplot(ranef(fm16.2Bmer))



# Test for (variance of) random slopes
#-------------------------------------#
# Auxiliary model with random slopes only. 
mAux <- update(fm16.1, random = ~ 0 + time|subject, data = armd)  

exactRLRT(m = mAux,      # Auxiliary model (the reduced model containing only the random effect to be tested)
          m0 = fm16.1,   # M16.1  (null)
          mA = fm16.2B)  # M16.2B (alternative)

# The simulated p-value is essentially 0
# indicating that null hypothesis can be rejected that the variance of random slopes is equal to 0.

# exactRLRT() function allows ONLY for Independent Random Effects.
# -->  we illustrate the use of the function for model M16.2B with a diagonal matrix D. (not for M16.2A)
# Because we consider a model with two variance components,
# i.e., random intercepts and random slopes, we need to specify all three arguments
# m, m0, and mA of the function exactRLRT()

?exactRLRT

