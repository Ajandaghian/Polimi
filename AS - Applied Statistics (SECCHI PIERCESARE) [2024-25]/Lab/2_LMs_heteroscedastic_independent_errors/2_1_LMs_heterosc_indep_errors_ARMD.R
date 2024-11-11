#__________________________ Applied Statistics 2024/2025 _________________________

#### 2.1.Linear Models with heteroscedastic and independent errors: ARMD Trial ####
#_________________________________________________________________________________#

# Topic:
# Linear Models with heteroscedastic and independent errors (theory: Chapter 7)
# in R: Chapters 8 and 9 for ARMD Trial dataset

# Set the directory!


#_________________________________________________________________________________________#
#### Fitting Linear Models with Heterogeneous Variance (& independent residual errors) ####
#_________________________________________________________________________________________#

library(nlme)  # for the models
library(nlmeU) # for the data
library(lattice)
library(corrplot)
library(plot.matrix)


# (Chapter 8.4)
# gls() function allows the inclusion of heteroscedasticity (and dependency, next lab)

# gls(model, 
#     data, 
#     subset,                           # optional
#     na.action,                        # optional
#     weights = varFunc(form=formula),  # Focus of this lab
#     control = glsControl()            # a list of control values to replace the default ones
                                        # e.g., for changing number of iterations, ...              
# )

# control
args(glsControl)
?glsControl


# weights
# NB. The default value of the weights argument is NULL --> LM with homoscedastic residual errors
# weights can be given directly as a one-sided formula

# varFunc Class (Chapter 8.2)
?varClasses

#_varClass___________parameters___________________Group
# varFixed()         value                        known weights
# varIdent()         value, form, fixed           <delta>-group
# varExp()           value, form, fixed           <delta>-group, <delta,mu>-group, <mu>-group
# varPower()         value, form, fixed           <delta>-group, <delta,mu>-group, <mu>-group
# varConstPower()    const, power, form, fixed    <delta>-group, <delta,mu>-group, <mu>-group


# what gls() returns
?glsObject   

# Extracting results directly from the object gls.fit = gls(...)
#_Model-fit component____________Syntax_________________________
# gls()-call                     (cl <- getCall(gls.fit))
# weights argument               cl$weights
# 95% CI for delta intervals     (gls.fit, which="var-cov")$varStruct
# Pearson residuals              resid(gls.fit, type="pearson")
# Var-cov structure              gls.fit$modelStruct
# Variance function              gls.fit$modelStruct$varStruct

?gls()  # see method

#_________________________________________________________________________________________#
#### ARMD Trial application (Chapter 9) : independent, heteroscedastic residual errors ####
#_________________________________________________________________________________________#
# heterogeneous variance, keeping assumption of independence

data(armd, package = "nlmeU")


# MODEL FORMULATION
# VISUAL_it = b_0t + b_1 * VISUAL0_i + b_2t * TREAT_i + eps_it
# eps_it~N(0,sigma_t^2), for t=1(4wks), 2(12wks), 3(24wks), 4(52wks)

# sqrt(Var(VISUAL_it)) = sigma_t where sigma_1=sigma_2=sigma_3=sigma_4
# so it means that it's the same as our previous lab: eps_it~N(0,sigma^2)

lm1.form <- visual ~ -1 + visual0 + time.f + treat.f:time.f   # model defined last time
fm6.1 <- gls(lm1.form, data = armd)
summary(fm6.1)

# with the command of last time we get the same result!
# fm6.1_lm <- lm(lm1.form, data = armd)
# summary(fm6.1_lm)

# Visualization of Variance-covariance matrix of Y (first 8 observations)
par(mar = c(4,4,4,4))
plot(diag(x=12.38^2, nrow=8, ncol=8), 
     col=colorRampPalette(c("blue", "white", "red")), 
     main='Variance-covariance matrix of Y')




##### varFixed() #####
?varFixed

###### model 9.0: lambda_i are known, i.e. lambda(v) ######

# we suppose that the variance is proportional to the time
# eps_it~N(0,sigma_t^2), for t=1(4wks),2(12wks),3(24wks),4(52wks)

#                                  {sigma*sqrt(4)   for t=1(4wks),
# sqrt(Var(VISUAL_it)) = sigma_t = {sigma*sqrt(12)  for t=2(12wks),
#                                  {sigma*sqrt(24)  for t=3(24wks),
#                                  {sigma*sqrt(52)  for t=4(52wks).
fm9.0 <- gls(visual ~ -1 + visual0 + time.f + treat.f:time.f, 
             method = 'REML',                   # default
             weights = varFixed(value = ~time), # Var.function; lambda_i(v_i) (v_i is known and observable)
             data = armd)
# NB. the variance covariate needs to be continuous: if we put time.f, it doesn't work!
# Try to change between 'REML' and 'ML' and see that results are different!
summary(fm9.0)

alpha = 0.05
intervals(fm9.0, which="coef", level = 1 - alpha)    # 1 - alpha/length(lm9.0$coefficients) if Bonferroni correction
intervals(fm6.1, which="coef", level = 1 - alpha)    # they are different!

intervals(fm9.0, which="var-cov", level = 1 - alpha)
intervals(fm6.1, which="var-cov", level = 1 - alpha) 

anova(fm6.1, fm9.0) # we can compare the models to see which one is better (lower AIC)

# Visualization of variance-covariance matrix of Y (8 observations, 2 patients)
par(mar = c(4,4,4,4))
plot(diag(x=c( 4 * 3.222976^2, # v_i * sigma^2
              12 * 3.222976^2, 
              24 * 3.222976^2, 
              52 * 3.222976^2), nrow=8, ncol=8),
     col=colorRampPalette(c("blue", "white", "red")),
     main='Variance-covariance matrix of Y - VarIdent() - Model 9.0')


# another example
# we suppose that the variance is proportional to 1/visual0 
fm9.0b <- gls(lm1.form, weights = varFixed(value = ~I(1/visual0)), data = armd)
summary(fm9.0b)

anova(fm6.1, fm9.0, fm9.0b)



##### varIdent() #####

###### model 9.1 ######
# <delta>-group i.e. lambda(delta,v)
# eps_it~N(0,sigma_t^2), for t=1(4wks),2(12wks),3(24wks),4(52wks)

#                                  {sigma*delta_1  for t=1(4wks),
# sqrt(Var(VISUAL_it)) = sigma_t = {sigma*delta_2  for t=2(12wks),
#                                  {sigma*delta_3  for t=3(24wks),
#                                  {sigma*delta_4  for t=4(52wks).

# delta_1 = 1
# delta_2 = sigma_2/sigma_1
# delta_3 = sigma_3/sigma_1
# delta_4 = sigma_4/sigma_1


lm1.form <- formula(visual ~ -1 + visual0 + time.f + treat.f:time.f)
fm9.1 <- gls(lm1.form,                  
         weights = varIdent(form = ~1|time.f), # Var.function; <delta>-group i.e. lambda(delta,v)
         data = armd)
# see ?varIdent for better understanding 'form':
# form ~ v, or ~ v | g, specifying a variance covariate v and, optionally, a grouping factor g for the coefficients.
# The variance covariate is ignored in this variance function.
# Defaults to ~ 1

summary(fm9.1)

fm9.1$modelStruct$varStruct             # delta1=1, delta_2 = sigma_2/sigma_1, ... 
# estimated values of the delta variance-function coefficients. 
# The estimates indicate an increasing variability of visual acuity measurements in time. 


# Visualization of variance-covariance matrix of Y (8 observations, 2 patients)
par(mar = c(4,4,4,4))
plot(diag(x=c(1.000000^2 * 8.244094^2, # delta_t^2 * sigma^2
              1.397600^2 * 8.244094^2, 
              1.664321^2 * 8.244094^2, 
              1.880852^2 * 8.244094^2), nrow=8, ncol=8),
     col=colorRampPalette(c("blue", "white", "red")),
     main='Variance-covariance matrix of Y - VarIdent() - Model 9.1')


# This is consistent with the results of exploratory analysis, e.g.:
bw1 <- bwplot(visual ~ time.f, data = armd)  
xlims <- c("4\nwks", "12\nwks", "24\nwks", "52\nwks")
update(bw1, xlim = xlims, pch = "|")

(intervals(fm9.1, which = "var-cov"))   # 95% CI for delta2,3,4 & sigma (see book 7.6.2)
# The standard deviation at week 4 is estimated to be equal to 8.24. 

# Comments:
# The 95% confidence intervals for the variance function coefficients slightly overlap, 
# but suggest timepoint-specific variances.


###### Comparison between the models ######

# To formally test the hypothesis that the variances are timepoint specific:
# --> Likelihood Ratio (LR) test statistic, based on the likelihood 
# for model fm6.1, which assumed homoscedasticity, and
# that for model fm9.1, which assumes heteroscedasticity. 

# Both models differ only by their variance structure and were fitted using the gls() function
# with the default estimation method (REML). 
# Moreover, model fm6.1 is nested within fm9.1, 
# because the former can be obtained from the latter by specifying 
# sigma_1^2 = sigma_2^2 = sigma_3^2 = sigma_4^2

# Thus, we can use the LR test to test the null hypothesis of homoscedasticity (fm6.1) 
anova(fm6.1, fm9.1)                     # fm6.1 nested in fm9.1
# pval < 0.0001 --> there is enough proof against the assumption H0 
# --> data provide evidence for heterogeneous variances of visual acuity measurements at != timepoints



##### varPower() ######

###### model 9.2 (REML estimation) ######
# <delta>-group i.e. lambda(delta,v)
# delta scalar, no strata 
# sigma_it = sigma * lambda_it 
#          = sigma * lambda(delta, TIME_it) 
#          = sigma * |TIME_it|^delta            
fm9.2 <- update(fm9.1,                    
         weights = varPower(form = ~time))           # <delta>-group 
# NB. continuous-time variable "time" rather than to the factor "time.f"
#     if you set time.f, there would be an error
?varPower # read form to better understand
fm9.2$modelStruct$varStruct                          # our delta is estimated

summary(fm9.2)

# Visualization of Variance-covariance matrix of Y (8 observations - 2 subjects)
par(mar = c(4,4,4,4))
plot(diag(x=c( 4^(2*0.2519332) * 5.974906^2, # TIME_it^{2*delta} * sigma^2
              12^(2*0.2519332) * 5.974906^2, 
              24^(2*0.2519332) * 5.974906^2, 
              52^(2*0.2519332) * 5.974906^2), nrow=8, ncol=8),
     col=colorRampPalette(c("blue", "white", "red")),
     main='Variance-covariance matrix of Y - VarPower() - Model 9.2')


###### Comparison between the models ######
anova(fm9.2, fm9.1)                        # fm9.2 nested in fm9.1: fm9.1 is more general
# Can a common-power variance function be used as a more parsimonious representation 
# of the variance structure of the data?
# we are testing 
# H0: sigma_1^2 = 4^delta  * sigma^2, 
#     sigma_2^2 = 12^delta * sigma^2, 
#     sigma_3^2 = 24^delta * sigma^2, 
#     sigma_4^2 = 52^delta * sigma^2
# vs
# H1: sigma_1^2 != 4^delta  * sigma^2, 
#     sigma_2^2 != 12^delta * sigma^2, 
#     sigma_3^2 != 24^delta * sigma^2, 
#     sigma_4^2 != 52^delta * sigma^2
#     (fm9.1)
# pval = 0.4 --> there is not enough proof/support against the assumption H0 
# --> fm9.2 is not statistically significantly worse than the fit of model fm9.1. 
# --> fm9.2, which specifies that the variance is a power function of the time (in weeks), 
#            offers an adequate description of the variance structure of the data.



###### model 9.3 (REML estimation) ######
# <delta>-group i.e. lambda(delta,v)
# delta = [delta_1, delta_2], strata = treatment group 

# sigma_it = sigma * lambda_it 
#          = sigma * lambda(delta, TIME_it) 
#          = {sigma * |TIME_it|^delta_1    ACTIVE
#            {sigma * |TIME_it|^delta_2    PLACEBO

fm9.3 <- update(fm9.1,                            
         weights = varPower(form = ~time|treat.f))   # <delta>-group
                                                     # strata=treat.f  
# NB. continuous-time variable "time" rather than to the factor "time.f"

fm9.3$modelStruct$varStruct
coef(fm9.3$modelStruct$varStruct)   # same, for getting coefficients

summary(fm9.3)
summary(fm9.3)$sigma


varWeights(fm9.3$modelStruct$varStruct)  # Variance weights: 1/^lambda_i
# inverse of the standard deviations corresponding to the variance function structure 
# for each observation in the dataset
                                        

# Visualization of Variance-covariance matrix of Y (8 random observations, 2 patients, one active, one placebo)
par(mar = c(4,4,4,4))
plot(diag(x=c(  4^(2*0.2532485) * 5.971599^2,  # ACTIVE - 4      # TIME_it^{2*delta_1} * sigma^2
               12^(2*0.2532485) * 5.971599^2,  # ACTIVE - 12
               24^(2*0.2532485) * 5.971599^2,  # ACTIVE - 24
               52^(2*0.2532485) * 5.971599^2,  # ACTIVE - 52
                4^(2*0.2511255) * 5.971599^2,  # PLACEBO - 4
               12^(2*0.2511255) * 5.971599^2,  # PLACEBO - 12
               24^(2*0.2511255) * 5.971599^2,  # PLACEBO - 24
               52^(2*0.2511255) * 5.971599^2), nrow=8, ncol=8), # PLACEBO - 52
     col=colorRampPalette(c("blue", "white", "red")),
     main='Variance-covariance matrix of Y - VarPower() - Model 9.3')


###### Comparison between the models ######
# Likelihood Ratio (LR) tests
anova(fm9.2, fm9.3)                        # fm9.2 nested in fm9.3
# we are testing 
# H0: delta_1=delta_2 (fm9.2)
# vs 
# H1: delta_1!=delta_2 (fm9.3)
# pval = 0.9 --> there is not enough proof/support against the assumption H0 (fm9.2)
# --> a common-power variance function of the TIME covariate can be used for both treatment groups.


###### model 9.4 (REML-based GLS) ######
# <delta,mu>-group i.e. lambda(delta,mu,v)
# delta = scalar, no strata 
# sigma_it = sigma * lambda_it 
#          = sigma * lambda(delta, mu_it) 
#          = sigma * |mu_it|^delta      where mu_it = b_0t + b1 * VISUAL0_i + b_2t * TREAT_i
#                                             predicted mean value of VISUAL_it

fm9.4 <- update(fm9.1, weights = varPower())                # <delta,mu>-group
# NB. equivalent to varPower(form=~fitted(.)), which means that the fitted values mu_i are used
# see ?varPower() for more details
fm9.4$modelStruct$varStruct

summary(fm9.4)

# Visualization of Variance-covariance matrix of Y (first 8 observations)
par(mar = c(4,4,4,4))
plot(diag(x=15.50761^2 * fitted(fm9.4)[1:8]^(-0.05885895*2), nrow=8, ncol=8),
     col=colorRampPalette(c("blue", "white", "red"))(10),
     main='Variance-covariance matrix of Y - VarPower() - Model 9.4')
# they are slightly different!
print(15.50761^2 * fitted(fm9.4)[1:8]^(-0.05885895*2))

###### Comparison between the models ######
# NOTICE. between the two best models fm9.4 and fm9.2, 
# we cannot use anova() (LRtest) since they are not nested!


###### model 9.5 (REML-based IRLS) ######
# <mu>-group i.e. lambda(mu,v)
# delta = 1, no strata 
# sigma_it = sigma * lambda_it 
#          = sigma * lambda(mu_it) 
#          = sigma * mu_it               
# NB. the scale parameter can be interpreted as a coefficient of variation, 
# constant for all timepoints.
fm9.5 <- update(fm9.1,                             
         weights = varPower(fixed = 1))              # <mu>-group
# NB. fixed=1 implies that, in the absence of any stratifying variable, 
# the power coefficient delta is fixed at 1.
# NB. as before, equivalent to varPower(form=~fitted(.)), 
# which means that the fitted values mu_i are used

fm9.5$modelStruct$varStruct
summary(fm9.5)

# Visualization of Variance-covariance matrix of Y (first 9 observations)
armd[1:9,] # 3 patients (3 "blocks": 11, 2222, 333)
par(mar = c(4,4,4,4))
plot(diag(x=0.2884781^2 * fitted(fm9.4)[1:9]^(1*2), nrow=9, ncol=9),
     col=colorRampPalette(c("blue", "white", "red"))(10),
     main='Variance-covariance matrix of Y - VarPower() - Model 9.5')
print(0.2884781^2 * fitted(fm9.4)[1:9]^(1*2))


###### Comparison between the models ######
anova(fm9.5, fm9.4)                        # fm9.5 nested in fm9.4 
# we are testing
# H0: if we assume a variance function in the form of a power function of the mean value, 
#     the power coefficient delta = 1 
#     (fm9.5)
# pval < 0.0001 --> reject the H0 (fm9.5) --> fm9.4 is more appropriate
# --> given that models fm9.4 and fm9.5 are both mean-variance models,
#     the inference on the implied variance structure, based on the result of the LR test, 
#     may need to be treated with caution (there are many unknowns)




##### Comparison between the models: AIC #####

# we check the AIC
AIC(fm9.1, fm9.2, fm9.3,                  
    fm9.4, fm9.5)                         # Smaller AIC corresponds to a better fit
# model fm9.2 offers a better fit to the data than model fm9.4
anova(fm9.2, fm9.4)                       # anova() command for Non-nested models

# NB. anova() is actually able to recognise whether two models are nested
anova(fm9.1, fm9.2, fm9.3,                  
      fm9.4, fm9.5)   

# --> we assess the fit of the best model (9.2) using residual plots


##### Residual analysis #####

library(lattice)

# Raw residuals 
plot(fm9.2, resid(., type = "response") ~ fitted(.))  # Raw residuals vs fitted
# We observe an asymmetric pattern, 
# with large positive (negative) residuals present mainly for small (large) fitted values.
# Our aim is to see no patterns at all in the residuals plot

plot(fm9.2, resid(., type = "response") ~ time)       # Raw residuals vs time 

bwplot(resid(fm9.2, type='response') ~ time.f,        # Raw residuals vs time.f
       pch = "|", data = armd)                        
# The box-and-whiskers plots clearly show an increasing variance of the residuals.



# Pearson residuals [ ^eps_i/sqrt(Var(y_i)) ]
# Pearson residuals are obtained from the raw residuals by dividing the latter by an
# estimate of the appropriate residual standard deviation, so they should be more homoscedastic

plot(fm9.2, resid(., type = "pearson" ) ~ fitted(.)) # Pearson vs fitted
# similarly to before, asymmetric pattern
plot(fm9.2, resid(., type = "pearson") ~ time)       # Pearson vs time 
bwplot( resid(fm9.2, type = "pearson") ~ time.f,     # Pearson vs time.f
        pch = "|", data = armd)
# thise plots illustrate the effect of scaling: 
# the variance of the residuals is virtually constant.


# Scale-location plots
# scale-location plots for the raw and Pearson residuals. 
# These are the scatterplots of the square-root transformation of the absolute
# value of the residuals versus fitted values. 
# --> The plots allow for detection of Patterns in the Residual Variance.

plot(fm9.2, sqrt(abs(resid(., type = "response"))) ~ fitted(.),
     type = c("p", "smooth"))
# there seems to be a dependence between the residual variance and the mean value. 
# However, this may be an artifact of the heteroscedasticity of the raw residuals

plot(fm9.2, sqrt(abs(resid(., type = "pearson"))) ~ fitted(.),
     type = c("p", "smooth"))
# does not indicate any clear trend in the residual variance







#### EXTRA material: NO NEED TO LEARN THIS CODE ####

# We select data for 188 subjects with all four post-randomization visual acuity measurements. 
# Scatterplot matrix of the Pearson residuals for model fm9.2 for all four measurement occasions.
# (complete cases only, n = 188; correlation coefficients above the diagonal)
# We also plot the 95% confidence ellipses.

residP <- resid(fm9.2, type="p")
dtAux <- subset(armd, select = c(subject, visual, time, time.f, treat.f)) 
require(reshape)
dtP <- data.frame(dtAux, residP)
dtPm <- melt(dtP,
             measure.var=c("residP"),
             id.var = c("subject","time.f"))
dtPc <- cast(subject ~ time.f ~ variable, data = dtPm) # array
dtPc <- data.frame(dtPc) 
names(dtPc) <- c("P4wks","P12wks","P24wks","P52wks")
head(dtPc)
range(dtPc, na.rm=TRUE)

library(ellipse)
my.upperPanel <-                           ## pairwise.complete.obs 
  function(x, y, subscripts, ...){
    panel.xyplot(x, y, type = "n", ...)      # no plot
    ic <- complete.cases(cbind(x,y))
    mn <- c(mean(x[ic]), mean(y[ic]))
    covx <- var(cbind(x,y), use="complete.obs")
    # print(covx)
    # ex <- ellipse(covx)
    corrx <- cov2cor(covx)
    corx <- round(corrx[1,2], 2)
    abs.corx <- abs(corx)
    # print(corx)
    cex.value <- 3
    ltext(0, 0, corx, cex = abs.corx * cex.value)
  }

my.lowerPanel <-                          ## pairwise.complete.obs 
  function(x,y,subscripts,...){
    panel.grid(h = -1, v = -1)
    covx <- var(cbind(x, y), use = "complete.obs")
    # print(covx)
    ex <- ellipse(covx)
    panel.xyplot(ex[ ,1], ex[ ,2], lty = 2, type = "l", ...)
    panel.xyplot(x, y, ...)
  }


mySuperPanel <- function(z, subscripts, panel.subscripts,...){
  panel.pairs(z, subscripts = subscripts,
              panel.subscripts = panel.subscripts,
              as.matrix=TRUE, 
              upper.panel = "my.upperPanel",
              lower.panel = "my.lowerPanel",
              prepanel.limits = function(z) return(c(-4,4))
  )
}


splom.form <- formula(~cbind(P4wks,P12wks,P24wks,P52wks))
splom.object <- splom(splom.form,
                      data=dtPc,             #### subset(armd240,miss.pat =="----"),
                      as.matrix=TRUE,        #### varnames = abbrev.names, 
                      xlab="",
                      superpanel = mySuperPanel
)

print(splom.object)
rm(my.upperPanel,mySuperPanel,splom.object)

# The scatterplots clearly show a violation of the assumption 
# of the independence of observations: residuals for different measurements are correlated. 
# The correlation coefficient decreases with the increasing distance between the timepoints. 

# NB. Some caution is needed in interpreting the strength of correlation, 
# because the estimated residuals are correlated even if the independence assumption holds


