#__________________________ Applied Statistics 2024/2025 _________________________

#### 1.1. Linear Models with homoscedastic* and independent errors: ARMD Trial ####
#_________________________________________________________________________________#
# *homoscedastic = homogeneous variance

# Topic: Linear Models with homoscedastic and independent errors (theory: Chapter 4)
# in R: Chapter 5 + Chapters 2.2 & 3.2 for ARMD Trial dataset + Chapter 6

# Summary:
# 1.1.1 Age-Related Macular Degeneration (ARMD) Trial: the Dataset (Chapter 2.2)
# 1.1.2 ARMD Trial: Data Visualization (Chapter 3.2)
# 1.1.3 ARMD Trial: Linear Model with Homogeneous Variance (Chapter 6) 
#                   [LM with independent, homoscedastic residual errors]

# SET THE DIRECTORY FIRST!

#________________________________________________________________________________________#
##### 1.1.1 Age-Related Macular Degeneration (ARMD) Trial: the Dataset (Chapter 2.2) #####
#________________________________________________________________________________________#

# A multi-center clinical trial comparing 
# an experimental treatment (interferon-alpha) and placebo in patients diagnosed with ARMD. 

# Patients with ARMD progressively lose vision. 
# Visual acuity of each of 240 patients participating in the trial was measured at baseline 
# and at (up to) 4 post-randomization timepoints (i.e., at 4, 12, 24, and 52 weeks). 

# Resulting data are an example of longitudinal data with observations grouped by subjects.

# Visual acuity was evaluated based on patient's ability to read lines of letters 
# on standardized vision charts. 
# The charts display lines of five letters of decreasing size, which the patient
# must read from top (largest letters) to bottom (smallest letters). 
# We will focus on the visual acuity defined as the total number of letters correctly read.

rm(list=ls())
library(settings)
reset(options)
graphics.off()

library(nlmeU)  # --> for the dataset
library(nlme)   # --> for models implementation
library(lattice)
library(corrplot)
library(plot.matrix)

data(armd.wide) # --> wide format
head(armd.wide) # lesion and line0 contain additional information, which we won't use
help(armd.wide)
dim(armd.wide)
str(armd.wide)

data(armd0)     # --> long (or longitudinal) format
head(armd0)     # new features: time.f, time, tp, and visual
help(armd0)
dim(armd0)
str(armd0)

data(armd)      # --> long (or longitudinal) format (subset of armd0, without baseline)
head(armd)
help(armd)
dim(armd)
str(armd)

# armd was built from armd0 as follows:
auxDt <- subset(armd0, time > 0)          # Post-baseline measures
dim(auxDt)                                # No. of rows & cols
levels(auxDt$time.f)                      # Levels of treat.f
armd <- droplevels(auxDt)                 # Drop unused levels   
levels(armd$time.f)                       # Baseline level dropped
armd <- within(armd, {                    # Contrasts assigned for dealing with ORDERED factors
  contrasts(time.f) <- contr.poly(4, scores = c(4, 12, 24, 52)) # we are changing scores
  # If you are interested in further details, see 1.5.1 in 1_5_Supplementary.R
})


#____________________________________________________________#
##### 1.1.2 ARMD Trial: Data Visualization (Chapter 3.2) #####
#____________________________________________________________#
# We are mainly interested in the effect of treatment on the visual acuity measurements.
# We plot measurements against time for several selected patients from both treatment groups. 

# Visual-acuity profiles for selected patients --> we visualize some of the trends

armd0.subset <- subset(armd0, as.numeric(subject) %in% seq(1, 240, 5)) # 1 every 5 patients

xy1 <- xyplot(visual ~ time | treat.f,   
              # visual and time are plotted against each other in separate
              # panels for different values of 'treat.f' factor.
              groups = subject,
              data = armd0.subset,
              type = "l", lty = 1)

update(xy1, xlab = "Time (in weeks)",
       ylab = "Visual acuity",
       grid = "h")

# Comments:
# - Decreasing trend in time, on average (patients progressively lose vision)
#   but patients have very different trends;
# - "Active" patients have on average lower values of the response;
# - There are patients for whom several measurements are missing;
# - Measurements adjacent in time seem to be well correlated, with the
#   correlation decreasing with an increasing distance in time;
# - Visual acuity at baseline seems to (partially) determine the overall level
#   of the post-randomization measurements.


# Inspecting missing-data patterns

# The factor miss.pat that indicates which of the four
# post-randomization measurements are missing for a particular patient
table(armd.wide$miss.pat)

# Comments:
# 188 patients for whom all four post-randomization visual acuity measurements were obtained.
# 6 patients for whom the four measurements were missing.
# 8 (= 4+1+2+1) patients with intermittent missing visual acuity measurements. 

# When modeling data with such patterns, extra care is needed 
# when specifying variance–covariance structures.


# Mean-Value profiles

# we compute the sample means 
# of visual acuity measurements for different visits and treatment groups (excluding na values)

# sample means across time and treatment
flst <- list(armd0$time.f, armd0$treat.f)                                     # "By" factors
(tN <-  tapply(armd0$visual, flst, FUN = function(x) length(x[!is.na(x)])))   # Counts

# NB: the function tapply() applies a selected function to each (nonempty) group of values 
# defined by a unique combination of the levels of one (or more) factors.

tMn <- tapply(armd0$visual, flst, FUN = mean)     # Sample means
tMd <- tapply(armd0$visual, flst, FUN = median)   # Sample medians
colnames(res  <- cbind(tN, tMn, tMd))             # Column names
nms1 <- rep(c("P", "A"), 3)
nms2 <- rep(c("n", "Mean", "Median"), rep(2, 3))
colnames(res) <- paste(nms1, nms2, sep = ":")     # New column names
res

# Comments:
# - on average, little difference in visual acuity between the two treatment groups at baseline. 
#   This is expected in a randomized study. 
# - During the course of the study, the mean visual acuity decreased with time 
#   in both settings (A and P)
# - Mean value is consistently higher in the placebo group, 
#   which suggests lack of effect of interferon-alpha.



# Box-plots for visual acuity for the five timepoints and the two treatments

bw1 <- bwplot(visual ~ time.f | treat.f, data = armd0)  # bwplot from package lattice
xlims <- c("Base", "4\nwks", "12\nwks", "24\nwks", "52\nwks")
update(bw1, xlim = xlims, pch = "|")

# Comments:
# - The decrease of the mean values in time is clearly seen for both treatment groups;
# - It is more pronounced for the active treatment arm;
# - As there was a slightly higher dropout in that arm, a possible
#   explanation could be that patients whose visual acuity improved dropped out of the study.


# Variance–covariance and correlation matrices for visual acuity

# measurements for complete cases only (n = 188)
visual.x <- subset(armd.wide, select = c(visual0:visual52))
(varx <- var(visual.x, use = "complete.obs"))            # Var-cov matrix
plot(varx, col=colorRampPalette(c("blue", "white", "red")))
# increase of the variance of visual acuity measurements obtained at later timepoints.
diag(varx)                                               # Var-cov diagonal elements

print(cor(visual.x, use = "complete.obs"), digits = 2)   # Corr matrix
cov2cor(varx)                                            # Corr matrix (alternative way) 
plot(cov2cor(varx), col=colorRampPalette(c("blue", "white", "red")))
# moderate-strong correlation which clearly decreases with the time gap

#______________________________________________________________________________#
##### 1.1.3 ARMD Trial: Linear Model with Homogeneous Variance (Chapter 6) #####
#______________________________________________________________________________#
data(armd, package = "nlmeU")

# NB: from now on, b stands for beta

# MODEL FORMULATION:
# VISUAL_it = b_0t + b_1 * VISUAL0_i + b_2t * TREAT_i + eps_it
# with eps_it~N(0,sigma^2)

# Comments:
# - VISUAL_it is the value of visual acuity measured for patient i (i = 1, . . . ,234) 
#     at time t (t = 1,2,3,4, corresponding to values of 4, 12, 24, and 52weeks, respectively);
# - b_0t is the timepoint-specific intercept;
# - b_1 is the baseline visual acuity effect;
# - VISUAL0_i is the baseline value of visual acuity;
# - b_2t is the timepoint-specific treatment effect.
# - TREAT_i is the treatment indicator (equal 1 for the active group and 0 otherwise);

# The model assumes a time-dependent treatment effect, 
# with the time variable being treated as a factor.

# MODEL FORMULATION in R:

lm1.form <- formula(visual ~ -1 + visual0 + time.f + treat.f:time.f ) 
# see 1_2_Formulation_Recap.R

# NB: To obtain timepoint-specific intercepts at 4,12,24 and 52 weeks, 
#    the overall intercept is removed from the model by specifying the -1 term.

# For further details related to design matrix and the contrasts for the (unordered) factors, 
# see 1.5.2 in 1_5_Supplementary_Materials.R

# we now fit the linear model
lm6.1 <- lm(lm1.form, data = armd)         # through lm()
lm6.1
summ <- summary(lm6.1)                     # Summary
summ

# NB: estimated coefficients for the time.f:treat.f interaction indicate
# a negative effect of the “Active” treatment.

# How to obtain information from our lm
lm6.1$coefficients       # the betas
lm6.1$residuals          # ^epsilon
lm6.1$fitted.values      # ^y
lm6.1$rank               # rank
lm6.1$df.residual        # degrees of freedom of the model (n-rank)
lm6.1$model              # the dataframe
rstandard(lm6.1)         # the standardized residuals

summ$sigma               # sigma (estimate of the residual standard deviation)

vcov(lm6.1)              # variance-covariance matrix
sqrt(diag(vcov(lm6.1)))  # standard errors of the coefficients


# Plot of the variance-covariance matrix of Y (first 8 observations) 
# --> it is a diagonal matrix with a value of 12.38^2
par(mar = c(4,4,4,4))
plot(diag(x=12.38^2, nrow=8, ncol=8), 
     col=colorRampPalette(c("blue", "white", "red")), 
     main='Variance-covariance matrix of Y')



# Confidence intervals for the coefficients
alpha = 0.05
confint(lm6.1,level = 1 - alpha)

# We can correct them through Bonferroni
confint(lm6.1,level = 1 - alpha/length(lm6.1$coefficients))

# We now try to use our model to fit new observations 
# and generate confidence and prediction intervals

# Step 1: we generate new observations - 3 patients (within the range of the observed values)
X.new <- data.frame(time.f = ordered(rep(c('4wks', '12wks', '24wks', '52wks' ),3)),
                    visual0 = c(rep(50,4), rep(60,4), rep(55,4)),
                    treat.f = factor(c(rep('Active',4), rep('Placebo',4), rep('Active',4)))
)
X.new
# Step 2: we build the intervals (for each single observation)
IC <-predict(lm6.1, X.new, interval="confidence", level=0.95)
IC
IP <-predict(lm6.1, X.new, interval="prediction", level=0.95)
IP
?predict.lm # look at the plot in the example below



anova(lm6.1)                               # ANOVA table
# Comment #1
# NB. p-value for the F-test for time.f:treat.f is significant (at the 5%)
# --> presence of a time-varying treatment effect. 
summary(lm6.1)
# the negative point estimates for the interaction coefficients indicate 
# a non-favorable treatment effect that increases over time. 

# Comment #2
# square of the value of the t-test statistic for visual0, 29.2132^2 = 853.43 
# does not equal the value of the F-test statistic, 14138.99
# why?

# the results of the t-tests provided by the summary() pertain to the MARGINAL testing.
# while anova() performs SEQUENTIAL approach (cfr book 4.7.1)

#       SEQUENTIAL APPROACH   MARGINAL APPROACH
# TEST  H0      H1            H0          H1
# X1    1       1+X1          1+X2+X3     1+X1+X2+X3
# X2    1+X1    1+X1+X2       1+X1+X3     1+X1+X2+X3
# X3    1+X1+X2 1+X1+X2+X3    1+X1+X2     1+X1+X2+X3

# t-test for visual0 assumes that visual0, time.f and time.f:treat.f, 
# are included in H0 (marginal)
# while the F-test in anova() assumes that no other terms 
# besides visual0 are present (sequential)

# Comment #3
# NB. anova() function can also be applied to more than one model-fit object
# providing the Likelihood Ratio tests for NESTED models
# as well as the values of AIC and BIC for each of the models

lm1.form2 <- formula(visual ~ -1 + visual0 + time.f)
lm6.1_2 <- lm(lm1.form2, data = armd)         # through lm()
summ2 <- summary(lm6.1_2)                     # Summary
summ2

anova(lm6.1_2, lm6.1)                         # fm6.1_2 nested in fm6.1
# look at Df (degrees of freedom added with model2)
# result of the test is statistically significant 
# --> model 2 (i.e. fm6.1) is better than model 1 (i.e. fm6.1_2)



# Hypothesis testing
# is it necessary to include in the model:...?
# 1. the variable treat.f;
# 2. the variable time.f;
# 3. the effect of the variable time.f on the intercept.
summary(lm6.1)

library(car)
# 1. the variable treat.f
linearHypothesis(lm6.1,
                 rbind(c(0,0,0,0,0,1,0,0,0),
                       c(0,0,0,0,0,0,1,0,0),
                       c(0,0,0,0,0,0,0,1,0),
                       c(0,0,0,0,0,0,0,0,1)),
                 c(0,0,0,0))

# 2. the variable time.f
linearHypothesis(lm6.1,
                 rbind(c(0,1,0,0,0,0,0,0,0),
                       c(0,0,1,0,0,0,0,0,0),
                       c(0,0,0,1,0,0,0,0,0),
                       c(0,0,0,0,1,0,0,0,0),
                       c(0,0,0,0,0,1,0,0,0),
                       c(0,0,0,0,0,0,1,0,0),
                       c(0,0,0,0,0,0,0,1,0),
                       c(0,0,0,0,0,0,0,0,1)),
                 c(0,0,0,0,0,0,0,0))
# 3. the effect of the variable time.f on the intercept
linearHypothesis(lm6.1,
                 rbind(c(0,1,0,0,0,0,0,0,0),
                       c(0,0,1,0,0,0,0,0,0),
                       c(0,0,0,1,0,0,0,0,0),
                       c(0,0,0,0,1,0,0,0,0)),
                 c(0,0,0,0))

# for further examples, see 1_6_Additional_exercises.R


#____________________________________________________#
# Diagnostics plots

# Scatterplot of raw residuals versus fitted values
plot(fitted(lm6.1), residuals(lm6.1)) # residuals() can be replaced by resid() or lm6.1$residuals      
abline(h = seq(-40, 40, by = 20), col = "grey")
abline(v = seq( 10, 80, by = 10), col = "grey")

# Comments: 
# The (vertical) width of the scatterplot clearly increases with increasing fitted values, 
# which implies a non-constant residual variance.

qqnorm(residuals(lm6.1))
qqline(residuals(lm6.1)) 

# Comments: the shape of the plot clearly deviates from a straight line. 
# Especially on the left tail.
# This may be an indication of a problem with the normality of the residuals. 
# However, it may also be the effect of ignored heteroscedasticity 
# and/or correlation of the visual acuity measurements.

shapiro.test(lm6.1$residuals) # problem confirmed by shapiro.test

# Automatic diagnostic plots (graphically)
x11()
par(mfrow = c(2,2))
plot(lm6.1)
# Third plot: checks whether residuals are spread equally along the ranges of predictors. 
# --> assumption of equal variance (homoscedasticity) is checked

# Residual analysis
plot(lm6.1$residuals) 
abline(h=0)
# homoscedastic?
# .. We remember that observations are not independent 
#    and that the variance of the visual measurements increases in time

# let's color the residuals relative to different patients
colori = rainbow(length(unique(armd$subject)))
num_sub = table(armd$subject)
colori2 = rep(colori, num_sub)
plot(lm6.1$residuals, col=colori2) 
abline(h=0)   
# --> different variability between patients

boxplot(lm6.1$residuals ~ armd$subject, col=colori,
        xlab='Subjects', ylab='Residuals', 
        main ='Distribution of residuals across patients')  # --> informative!

# let's color the residuals relative to different time instants
set.seed(1)
colori =rainbow(4)
colori2 = colori[armd$tp] # associate to each one of the 4 time instants a color
plot(lm6.1$residuals, col=colori2, ylab='residuals')
abline(h=0)
legend(650, -25, legend=c("time 4wks", "time 12wks", "time 24wks", "time 52wks"),
       col=colori, lty=1, cex=0.7) # --> not very informative!

# Comments: we observe that red points are the closest to 0, purple ones are the furthest
# We expect the residuals to be heterogeneous across different time instants observations

boxplot(lm6.1$residuals ~ armd$time.f, col=colori,
        xlab='Time.f', ylab='Residuals') # --> informative
# -> the variance of the observations increases in time

# The model does not take into account the correlation 
# between the visual acuity observations obtained from the same subject. 
# It also does not take into account the heterogeneous variability
# present at different time points. 
# Thus, it should not be used as a basis for inference.


#____________________________________________________#
##### The model can be also fitted through gls() #####
#____________________________________________________#
require(nlme)              # Attach nlme package
fm6.1 <- gls(lm1.form, data = armd)
summary(fm6.1)
intervals(fm6.1)           # 95% CI for beta, sigmas

plot(predict(fm6.1), residuals(fm6.1))   # Same as before
qqnorm(residuals(fm6.1))                 # Same as before
qqline(residuals(fm6.1))

# see 1_3_lm()_VS_glm()_Recap.R
# for the main differences between the two functions



