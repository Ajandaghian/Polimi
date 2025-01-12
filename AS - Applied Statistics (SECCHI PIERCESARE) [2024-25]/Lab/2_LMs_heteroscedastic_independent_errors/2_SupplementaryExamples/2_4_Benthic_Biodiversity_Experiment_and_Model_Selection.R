#__________________________ Applied Statistics 2024/2025 _________________________

#### 2.4.Benthic Biodiversity Experiment: 
####     more on LMs with heteroscedastic and independent errors and the Model Selection Procedure 
#_________________________________________________________________________________#

# The data relies on both published data (Ieno et al., 2006) and unpublished data (University of Aberdeen). 
# The experiment examines the effect of macrofauna density (H. diversicolor, Polychaeta), 
# and habitat heterogeneity on sediment nutrient release.

# At the start of the experiment, each container (n=108) is filled with homogenised sediment 
# from mudflats on the Ythan estuary (Scotland, UK). 
# The macrofaunal biomass (H. diversicolor) was fixed across levels 0, 0.5, 1, 1.5, and 2 g. 
# The response variable is the concentration of a particular nutrient (PO3, NO3, NH4).
# --> 36 observations per 3 nutrients, 18 enriched (Algae), and 18 non-enriched (NoAlgae) 
# --> tot 108

# Source: Mixed Effects Models and Extensions in Ecology with R (2009) Zuur, Ieno, Walker, Saveliev and Smith. Springer

Biodiv = read.table('Biodiversity.csv', header = TRUE, dec=".", sep=',')

str(Biodiv)
Biodiv$fTreatment = as.factor(Biodiv$Treatment)
Biodiv$fNutrient = as.factor(Biodiv$Nutrient)

# Due to the nature of the variables, we expect massive differences in variation
# in concentrations per nutrient and enrichment combination.
# Boxplot for each nutrient–enrichment combination;
boxplot(Concentration ~ fTreatment * fNutrient, data = Biodiv)
# The samples enriched with algae and with NH4, 
# have higher concentrations and show more variation

# Let's fit our first "classical" model
# --> Linear regression analysis, using biomass, enrichment and nutrient, with
# all the two-way interactions, and the three-way interaction as explanatory variables

mod0 = lm(Concentration ~ Biomass * fTreatment * fNutrient, data=Biodiv)
summary(mod0)

# Automatic diagnostic plots (graphically)
par(mfrow = c(2,2))
plot(mod0)

shapiro.test(mod0$residuals)

# Clearly shows serious violation of homogeneity
# Biological knowledge suggests that treatment and nutrient levels, possibly both, 
# may be driving the heterogeneity. 

# Let's visualize it

# Residuals & Nutrient
set.seed(1)
colori = rainbow(length(levels(Biodiv$fNutrient)))
colori2 = colori[Biodiv$fNutrient] 
plot(mod0$residuals, col=colori2, ylab='residuals')
abline(h=0)

boxplot(mod0$residuals ~ Biodiv$fNutrient, col=colori, xlab='Nutrient', ylab='Residuals') 

# Residuals & Treatment
colori =rainbow(length(levels(Biodiv$fTreatment)))
colori2 = colori[Biodiv$fTreatment] 
plot(mod0$residuals, col=colori2, ylab='residuals')
abline(h=0)

boxplot(mod0$residuals ~ Biodiv$fTreatment, col=colori, xlab='Treatment', ylab='Residuals')


plot(Biodiv$Biomass, Biodiv$Concentration)
# No clear increase or decrease in spread. 
# --> the potential variance covariates are nutrient and/or enrichment.

#_______________________________________________________________________________
# Let us fit now some models by introducing a variance structure
library(nlme)
# linear model without any variance covariates
M0 <- gls(Concentration~Biomass*fTreatment*fNutrient, data = Biodiv)
plot(M0, col = 1)

M1A <- gls(Concentration~Biomass*fTreatment*fNutrient,
         weights=varIdent(form=~1|fTreatment*fNutrient),
         data = Biodiv) # just one variance covariate per nutrient-enrichment combination

M1B<-gls(Concentration~Biomass*fTreatment*fNutrient,
         weights=varIdent(form=~1|fNutrient),
         data=Biodiv)   # Nutrient covariates

M1C<-gls(Concentration~Biomass*fTreatment*fNutrient,
         weights=varIdent(form=~1|fTreatment),
         data = Biodiv) # Enrichment covariates


anova(M0, M1A, M1B, M1C)
# Looks like varIdent model with nutrient-enrichment combination (M1A) is best (see AIC).

summary(M1A)
plot(M1A, col = 1) # this command plots the standardized residuals versus fitted values
# No sign of heterogeneity;
# We have found the optimal residual variance structure using REML estimation.

#_______________________________________________________________________________
#### Which explanatory variables are significant? ####

# 3 tools for variables selection:
# 1) t-statistic -> summary(model), model fitted through 'REML'
#                   NB. t-statistic is should not be used to assess the significance 
#                       of a categorical variable with more than two levels.
# 2) F-statistic -> anova(model), model fitted through 'REML'; anova does SEQUENTIAL TESTING.
#                   NB. this is  useful for testing the significance 
#                       of the highest interaction term, but not for the other terms in the model. 
#                       It is also less useful if you only have main terms 
# 3) Likelihood ratio test -> we need to specify a full model and a nested model. 
#                             NB. Both models need 'ML' estimation 
#                                 (and the "optimal" variance function we found) 

summary(M1A) # from the summary we get t-statistics
anova(M1A)   # We are performing SEQUENTIAL TESTING
# This means that the p-values will change if the order of the main terms 
# or the order of the two-way interactions is changed
# e.g. model: beta0 + beta1 * x1 + beta2 * x2 + eps
#             H0           H1         model
# beta0:   beta0 = 0   beta0 != 0     beta0
# beta1:   beta1 = 0   beta1 != 0     beta0 + beta1 * x1
# beta2:   beta2 = 0   beta2 != 0     beta0 + beta1 * x1 + beta2 * x2 

# only the last term that is of real interest as it shows the significance
# of the three-way interaction term (order of this term cannot be changed)
# both functions show that the three-way interaction is not significant:
# we can drop the three-way term and refit the model.

# Let us check also through the likelihood ratio test method (with 'ML')
# if we can drop the three way interaction:

M2A1 <- gls(Concentration ~ Biomass + fTreatment + fNutrient +
            Biomass:fTreatment + Biomass:fNutrient + fTreatment:fNutrient +
            Biomass:fTreatment:fNutrient,
          weights = varIdent(form =~ 1 | fTreatment * fNutrient),
          method = "ML", data = Biodiv)

M2A2 <- gls(Concentration ~ Biomass + fTreatment + Nutrient +
            Biomass:fTreatment + Biomass:fNutrient + fTreatment:fNutrient,
          weights=varIdent(form =~ 1 | fTreatment * fNutrient),
          method="ML", data = Biodiv)

anova(M2A1, M2A2)

# Also the anova command indicates that the three-way interaction can be dropped. 

# What about the two-way interactions?
# 1) t-statistic
# we cannot use t-statistic because nutrient has three levels. 
# One level will be used as baseline, and the p-values from the t-statistic will
# only tell us whether the second and third nutrients are different from the baseline nutrient.

# 2) F-statistic
# with anova, we cannot assess the significance of the Biomass*Treatment term, 
# and the biomass*Nutrient term, due to the order how they were put in. 
# [We could apply 3 models, ensure each time that a different two-way term is the last, 
# and deselect the least significant two-way interaction]

# 3) Likelihood-ratio test
# We can compare nested models.
# NB. Both models need ML estimation. 

# -> we proceed with method 3)

#_______________________________________________________________________________
##### Round 1 of the Backwards Selection #####
fFull <- formula(Concentration ~ Biomass + fTreatment + fNutrient +
                                 Biomass:fTreatment + Biomass:fNutrient + fTreatment:fNutrient)
vfOptim <- varIdent(form =~ 1 | fTreatment*fNutrient)

M3.Full<-gls(fFull,
             weights = vfOptim,
             method = "ML", data = Biodiv)  # with Backward Selection, remember to use "ML" estimation

# Drop Biomass:fTreatment?
M3.Drop1<-update(M3.Full,.~.-Biomass:fTreatment)
anova(M3.Full, M3.Drop1) # p-val > 0.05

# Drop Biomass:fNutrient?
M3.Drop2<-update(M3.Full,.~.-Biomass:fNutrient)
anova(M3.Full, M3.Drop2) # p-val < 0.05

# Drop fTreatment:fNutrient?
M3.Drop3<-update(M3.Full,.~.-fTreatment:fNutrient)
anova(M3.Full,M3.Drop3)  # p-val < 0.05

# --> we drop the term Biomass:fTreatment

#_______________________________________________________________________________
##### Round 2 of the Backwards Selection #####
# New full model
M4.Full<-gls(Concentration ~ Biomass + fTreatment + fNutrient +
                             Biomass:fNutrient + fTreatment:fNutrient,
             weights= vfOptim,
             method="ML", data=Biodiv)
# From this model, we can drop two of the two-way interaction terms. 
# No main terms can be dropped yet!

# Drop Biomass:fNutrient?
M4.Drop1 <- update(M4.Full, .~. -Biomass:fNutrient )
anova(M4.Full, M4.Drop1)     # A p-value of 0.04 for the Biomass:fNutrient interaction 
                             # is not impressive, especially not with a series 
                             # of hypothesis tests. So, we drop it 

# Drop fTreatment:fNutrient?
M4.Drop2 <- update(M4.Full, .~. -fTreatment:fNutrient )
anova(M4.Full, M4.Drop2)     # p-val < 0.05

# --> we drop the term Biomass:fNutrient

#_______________________________________________________________________________
##### Round 3 of the Backwards Selection #####
# New full model
M5.Full<-gls(Concentration ~ Biomass + fTreatment + fNutrient + fTreatment:fNutrient,
             weights=vfOptim,
             method="ML", data=Biodiv)

# The rule is that if an interaction term is included, 
# then all the associated main terms should be included as well, 
# and are not a candidate for dropping. 
# However, if you have the main terms A, B, C, and the interaction A:B, 
# then the two terms that can be potentially dropped are A:B and also C!

# Drop fTreatment:fNutrient?
M5.Drop1 <- update(M5.Full, .~. -fTreatment:fNutrient)
anova(M5.Full, M5.Drop1)     # p-val < 0.05

# Drop Biomass (main term)?
M5.Drop2 <- update(M5.Full, .~. -Biomass)
anova(M5.Full, M5.Drop2)     # p-val > 0.05

# --> we drop the term Biomass

#_______________________________________________________________________________
##### Round 4 of the Backwards Selection #####
# New full model
M6.Full<-gls(Concentration ~ fTreatment + fNutrient + fTreatment:fNutrient,
             weights=vfOptim,
             method="ML", data=Biodiv)

# Drop fTreatment:fNutrient?
M6.Drop1 <- update(M6.Full, .~. -fTreatment:fNutrient)
anova(M6.Full, M6.Drop1)   # p-val < 0.05

# --> we do not drop anything

#_______________________________________________________________________________
##### Final model #####
# Once we understood which features to drop, we refit the model with 'REML' estimation
MFinal<-gls(Concentration ~ fTreatment + fNutrient + fTreatment:fNutrient,
            weights=vfOptim,
            method="REML", data=Biodiv)

# Can homogeneity be safely assumed?

E <- resid(MFinal, type="pearson")
Fit = fitted(MFinal)

# Note on resid() function:
# If "response", as by default, the “raw” residuals (observed - fitted) are used; 
# else, if "pearson", the standardized residuals (raw residuals divided by 
#                                                 the corresponding standard errors) are used; 
# else, if "normalized", the normalized residuals (standardized residuals pre-multiplied 
#                                                  by the inverse square-root factor 
#                                                  the estimated error correlation matrix) are used. 


op <- par(mfrow=c(1,2))
plot(x=Fit, y=E, xlab="Fitted values", ylab="Std residuals", main="Std residuals versus fitted values")
# identify(Fit, E)   # identify the observation with the large residual (observation 26)
hist(E, nclass=15)   # a histogram of the residuals (denoted by E) for the optimal GLS model

# a comparison between the 'classical' model and the obtained one
plot(M0, col = 1)
plot(MFinal, col = 1)

library(lattice)
# Pearson residuals on M0
plot(M0, resid(., type = "pearson") ~ fitted(.))  # Raw residuals vs fitted
bwplot(resid(M0, type='pearson') ~ fNutrient,     # Raw residuals vs fNutrient
       pch = "|", data = Biodiv)                        
bwplot(resid(M0, type='pearson') ~ fTreatment,    # Raw residuals vs fTreatment
       pch = "|", data = Biodiv)  


# Pearson residuals on MFinal
plot(MFinal, resid(., type = "pearson") ~ fitted(.))  # Raw residuals vs fitted
# the variance of the residuals is virtually constant. 
bwplot(resid(MFinal, type='pearson') ~ fNutrient,     # Raw residuals vs fNutrient
       pch = "|", data = Biodiv)                        
bwplot(resid(MFinal, type='pearson') ~ fTreatment,    # Raw residuals vs fTreatment
       pch = "|", data = Biodiv)  

# scale-location plots: M0 VS MFinal
plot(M0, sqrt(abs(resid(., type = "pearson"))) ~ fitted(.),
     type = c("p", "smooth"))

plot(MFinal, sqrt(abs(resid(., type = "pearson"))) ~ fitted(.),
     type = c("p", "smooth"))
# does not indicate any trend in the residual variance


# Pearson residuals for MFinal:
# Residuals & Nutrient
set.seed(1)
colori =rainbow(length(levels(Biodiv$fNutrient)))
colori2 = colori[Biodiv$fNutrient] 
plot(E, col=colori2, ylab='residuals')
abline(h=0)

boxplot(E ~ Biodiv$fNutrient, col=colori, xlab='Nutrient', ylab='Residuals') 

# Residuals & Treatment
colori =rainbow(length(levels(Biodiv$fTreatment)))
colori2 = colori[Biodiv$fTreatment] 
plot(E, col=colori2, ylab='residuals')
abline(h=0)

boxplot(E ~ Biodiv$fTreatment, col=colori, xlab='Treatment', ylab='Residuals')


# Can normality be safely assumed?
qqnorm(residuals(MFinal))
qqline(residuals(MFinal)) 
shapiro.test(MFinal$residuals)


# Let's have a look at the summary
summary(MFinal) 
# Note that the combination enrichment with algae and NH4 has the largest variance, 
# namely (8.43 * 0.819)^2.
# All terms are significantly different from 0 at the 5% level. 

# To understand what the model telling us,
# it can be helpful to graph the fit of the model.
boxplot(predict(MFinal) ~ fTreatment:fNutrient, data = Biodiv) 
# This only works because all the explanatory variables are nominal.
# The observations exposed to algae treatment and NH4 enrichment have the highest values. 
# This explains why the interaction term is significant. 


# Homework:
# Try to remove row 26 from the data, or add subset = –26 to each gls command
# to see what happens


